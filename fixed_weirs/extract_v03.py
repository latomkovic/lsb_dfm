# v03: use fixed_weirs-v01.shp, to bring in some manual features
#      also means that we no longer need to filter things out
# v02: pull max elevation within 5m -- at some places the levee line does not
# quite follow the crest of the levee.  
# 
# v01: separate out by 'Class'

import os
import matplotlib.pyplot as plt
from matplotlib import cm

from shapely import geometry
from scipy import ndimage

from stompy.spatial import (field,wkb2shp,linestring_utils)
import stompy.plot.cmap as scm
from stompy.grid import unstructured_grid
from stompy.model.delft import io

data_root='/opt/data'

## 

dem=field.GdalGrid(os.path.join(data_root,'bathy_interp/master2017/tiles_2m_20171024/merged_2m.tif'))

## 

# More details on the provenance of this on the Google drive in
# Modeling/Lower South Bay/shoreline elevs
inv_fn="fixed_weirs-v01.shp"
inv=wkb2shp.shp2geom(inv_fn)

## 

# For the moment, focus on just the parts of the inventory which overlap with the DEM.
bounds_xyxy=np.array(dem.extents)[ [0,2,1,3] ]
box=geometry.box(*bounds_xyxy)

sel=[ box.intersects(geo) for geo in inv['geom']]

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

pli_data=io.grid_to_pli_data(g,node_fields=['elev_m','sill_left','sill_right'],
                             labeler=lambda i: "L%04d"%(total_count+i))
io.write_pli('fixed_weirs-v01.pli',pli_data)

## 

