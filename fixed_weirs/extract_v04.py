# v04: Add code to process water control structures.
# v03: use linear_features-v02.shp, to bring in some manual features,
#      mostly these are the features we want, although gate-like features
#      are also in here.
# v02: pull max elevation within 5m -- at some places the levee line does not
# quite follow the crest of the levee.  
# v01: separate out by 'Class'

import os
import matplotlib.pyplot as plt
from matplotlib import cm

import numpy as np

from shapely import geometry
from scipy import ndimage

from stompy.spatial import (field,wkb2shp,linestring_utils)
from stompy import utils
import stompy.plot.cmap as scm
from stompy.grid import unstructured_grid
import stompy.model.delft.io as dio

data_root='..'
out_dir='out'

## 

os.path.exists(out_dir) or os.makedirs(out_dir)
## 

dem=field.GdalGrid(os.path.join(data_root,'inputs-static/merged_2m.tif'))

## 

# More details on the provenance of this on the Google drive in
# Modeling/Lower South Bay/shoreline elevs
inv_fn="linear_features-v02.shp"
inv=wkb2shp.shp2geom(inv_fn)

## 

# For the moment, focus on just the parts of the inventory which overlap with the DEM.
bounds_xyxy=np.array(dem.extents)[ [0,2,1,3] ]
box=geometry.box(*bounds_xyxy)

# There are some water control structures which came in from the original
# shoreline inventory.  Only structures which have been named are omitted here -
# these are the ones which will be handled specifically as gates, while other
# structures will be left as closed
sel=[ box.intersects(geo) and ( (klass!='Water Control Structure') or (model_name==''))
      for model_name,klass,geo in zip(inv['model_name'],inv['Class'],inv['geom'])]

# trims by a factor of 8
sel=np.array(sel)

## 

# New in v02:
# Apply a max filter across the DEM.

# for a 5m radius and 2m pixels, not much hope in really resolving
# a disc, but here goes
footprint=np.array( [[0,1,1,1,0],
                     [1,1,1,1,1],
                     [1,1,1,1,1],
                     [1,1,1,1,1],
                     [0,1,1,1,0]] )

## 

# The real deal - update dem in place (in RAM, not on disk)
dem.F=ndimage.maximum_filter(dem.F,footprint=footprint)

## 
res=field.ConstantField(5.0) # target output linear resolution

# count total features so that files can be concatenated without issue.
total_count=0

# no longer split by class, all filtering already complete in the input
# shapefile

g=unstructured_grid.UnstructuredGrid(extra_node_fields=[ ('elev_m','f8')],
                                     extra_edge_fields=[ ('mean_elev_m','f8')] )

for ix,sel_i in enumerate(np.nonzero(sel)[0]):
    geom=inv['geom'][sel_i]
    coords=np.array(geom)
    # not 100% sure about whether it's necessary to test for closed_ring
    new_coords=linestring_utils.upsample_linearring(coords,res,closed_ring=0)

    nodes=[g.add_or_find_node(x=xy,tolerance=0.0,
                              elev_m=np.nan)
           for xy in new_coords]

    for a,b in zip(nodes[:-1],nodes[1:]):
        j=g.add_edge(nodes=[a,b],mean_elev_m=inv['Z_Mean'][sel_i])

# pull out point elevations at the nodes:
g.nodes['elev_m'] = dem( g.nodes['x'] )

# drastic, but go ahead and delete any nodes which failed to get an elevation

missing=np.isnan( g.nodes['elev_m'] )

for n in np.nonzero(missing)[0]:
    g.delete_node_cascade(n)

## 

# to replicate the 5 fields, where last two are just 10.0
# not entirely sure what these /should/ be, but this is what
# I've seen in previous input files.
g.add_node_field('sill_left',10*np.ones_like(g.nodes['elev_m']))
g.add_node_field('sill_right',10*np.ones_like(g.nodes['elev_m']))

## 

pli_data=dio.grid_to_pli_data(g,node_fields=['elev_m','sill_left','sill_right'],
                             labeler=lambda i: "L%04d"%(total_count+i))
dio.write_pli(os.path.join(out_dir,'fixed_weirs-v02.pli'),pli_data)

## 


def write_node_shp(self,shpname,extra_fields=[]):
    """ Write a shapefile with each node.  Fields will attempt to mirror
    self.nodes.dtype

    extra_fields: goal is similar to write_cells_shp and write_edges_shp, 
    but not yet supported.
    """
    assert len(extra_fields)==0 # not yet supported!
    
    # zero-based index of node (why does write_edge_shp create 1-based ids?)
    base_dtype = [('node_id',np.int32)]

    node_geoms=[geometry.Point( self.nodes['x'][i] )
                for i in self.valid_node_iter() ]

    node_data=self.nodes[~self.nodes['deleted']].copy()

    # don't need to write all of the original fields out:
    node_data=utils.recarray_del_fields(node_data,['x','deleted'])
    

    wkb2shp.wkb2shp(shpname,input_wkbs=node_geoms,fields=node_data,
                    overwrite=True)


nodes_fn=os.path.join(out_dir,inv_fn.replace(".shp","-nodes.shp"))
write_node_shp(g,nodes_fn)

##

# Gates!

# Only structures which have been named are included here -
# other structures will be left as closed, and handled above.
sel=[ (klass=='Water Control Structure') and (model_name!='')
      for model_name,klass,geo in zip(inv['model_name'],inv['Class'],inv['geom'])]
sel=np.array(sel) 
idxs=np.nonzero(sel)[0]


# 6 known gates:
# a5_guad, a7_alviso, a9_alviso, a14_coyote, a8_alviso, a3w_guad
# For starters, all of these get the same, constant in time opening
# the resulting ini file is the StructureFile referred to in the mdu
with open(os.path.join(out_dir,'gates-v04.ini'),'wt') as fp:
    for idx in idxs:
        pli_base_fn="gate-%s.pli"%inv['model_name'][idx]
        pli_fn=os.path.join(out_dir,pli_base_fn)
        dio.write_pli(pli_fn, [ (inv['model_name'][idx], # label
                                 np.array(inv['geom'][idx]) ) ] )

        # these parameters are a good starting point for alviso_a8
        # may be too wide and/or too shallow for others.

        # not entirely sure of how these geometries are interpreted
        # Okay - after much trial and error, seems like
        # sill_level is the elevation of the fixed, weir-like structure below
        # a moving gate.
        # lower_edge_level is the elevation of a gate which may be flush with
        # sill level (flow is allowed only by horizontal opening), below the
        # sill level (functionally the same as flush), or above the sill_level
        # (flow is allowed through a vertical gap between sill and the gate)
        # opening_width: the horizontal opening, i.e. a gate that slides horizontally
        # and partially blocks the flow.  NOTE!!! this appears to work only in
        # integer number of flux faces.  i.e. opening_width cannot be used to
        # create a 1m wide sluice, unless you have 1m wide grid cells.
        # door_height: should probably just be very large.  In theory this would
        # allow flow over the gate, but some comments in the code suggest that
        # is not implemented.
        # sill_width: an upper bound on opening width, defaults to the length of
        # the flux faces

        name=inv['model_name'][idx]
        
        fp.write("[structure]\n")
        fp.write("type                         = gate\n")
        fp.write("id                           = %s\n"%name)
        fp.write("polylinefile                 = %s\n"%pli_base_fn)
        fp.write("door_height                  = 15\n")
        
        if name=='a8_alviso':
            fp.write("lower_edge_level             = 1\n") # This can be a tim file
            fp.write("opening_width                = 25\n") # This can be a tim file
            fp.write("sill_level                   = 1\n")
        else:
            # More like a culvert
            fp.write("lower_edge_level             = -0.95\n") # This can be a tim file
            fp.write("opening_width                = 0\n") # This can be a tim file
            fp.write("sill_level                   = -1.0\n")
            
        # the GUI seems to require this to be present.
        fp.write("horizontal_opening_direction = symmetric\n")
        fp.write("\n")
