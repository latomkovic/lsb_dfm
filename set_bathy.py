"""
The incoming grid has no bathymetry.

This script sets bathymetry on the nodes in two ways:
1. A first pass interpolates from the sfb_dfm_v2 grid, so most nodes
   will get the same value since most of the nodes are unchanged from
   that grid.
2. A second pass evaluates the 2m LSB bathymetry on edges, then uses a 
   "greedy" strategy to move those values to nodes, assuming bed level type
   of 4 (4, right?)

This script relies on file locations on hpc, as of 2017-08-30, for 
the 2m LSB bathy (merged_2m_path below).

The original grid bathy is taken from the sfb_dfm submodule.

This module defines the method set_lsb_bathy(grid), but it is up to
a caller to invoke this (i.e. lsb_dfm.py)
"""

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from stompy import utils
from stompy.grid import unstructured_grid
from stompy.model.delft import dfm_grid
from stompy.spatial import field
from stompy.plot import plot_utils
from stompy.grid import depth_connectivity

import logging
log=logging.getLogger('set_bathy')

##

merged_2m_path='inputs-static/merged_2m.tif'

if not os.path.exists(merged_2m_path):
    raise Exception("Copy or symlink merged_2m.tif to %s"%merged_2m_path)


def set_lsb_bathy(grid):
    if 0:
        plt.figure(10).clf()
        fig,ax=plt.subplots(num=10)
        ecoll=g.plot_edges(lw=0.3,color='k')
        ncoll=g.plot_nodes(values=g.nodes['depth']) # ds.NetNode_z
        plot_utils.cbar(ncoll)

    # First, load in the original sfb_dfm grid to get bathymetry
    log.info("Loading SFB DFM v2 grid")
    sfb_dfm_grid=xr.open_dataset('sfb_dfm/sfei_v19_net.nc')

    sfb_X=np.c_[ sfb_dfm_grid.NetNode_x.values,
                 sfb_dfm_grid.NetNode_y.values]
    sfb_dfm_field=field.XYZField(X=sfb_X,
                                 F=sfb_dfm_grid.NetNode_z.values)

    lsb_X=grid.nodes['x']

    # This will set all points within the convex hull of the original
    # grid.  These elevations are used in the output (a) outside the
    # LSB merged_2m DEM, and (b) to prioritize which nodes to use when
    # moving depths from edges to nodes.
    # Could argue that it would be better to pull point elevations here
    # from the DEM where they overlap.  Probably makes very little difference
    lsb_z=sfb_dfm_field(lsb_X) # still some nans at this point

    #
    
    log.info("Loading LSB dem from %s"%merged_2m_path)
    dem=field.GdalGrid(merged_2m_path)

    #

    # Use the 2m DEM to find an effective minimum depth for each edge covered
    # by the DEM.
    edge_min_depths=depth_connectivity.edge_connection_depth(grid,dem,centers='lowest')
    # move those edge depths to node depths
    node_depths=depth_connectivity.greedy_edgemin_to_node(grid,lsb_z,edge_min_depths)

    # Still have some nodes with nan depth, first try filling in with the DEM.
    missing=np.isnan(node_depths)
    node_depths[missing]=dem( lsb_X[missing,:] )

    # And wrap it up with a more forgiving interpolation from the original sfb_dfm_v2
    # grid (about 12 points on the convex hull)
    still_missing=np.isnan(node_depths)
    node_depths[still_missing]=sfb_dfm_field.interpolate( lsb_X[still_missing,:],
                                                          'nearest' )
    assert np.isnan(node_depths).sum()==0
    
    # Update the grid
    grid.nodes['depth']=node_depths

    if 0: # caller is going to deal with I/O
        out_file='lsb_v99_bathy_net.nc'
        os.path.exists(out_file) and os.unlink(out_file)
        dfm_grid.write_dfm(grid,out_file)

    if 0: # plot for development.
        plt.figure(10).clf()

        fig,ax=plt.subplots(num=10)

        edge_mins=grid.nodes['depth'][grid.edges['nodes']].min(axis=1)

        ecoll=grid.plot_edges(lw=1.5,values=edge_mins)
        ncoll=grid.plot_nodes(values=grid.nodes['depth'])
        plt.setp([ecoll,ncoll],clim=[-3,3])
        plot_utils.cbar(ncoll,extras=[ecoll])

    # Modified in place, but return just in case
    return grid # QED.

##

if 0:
    g=dfm_grid.DFMGrid('lsb_v99_net.nc')
    g2=set_lsb_bathy(g)
