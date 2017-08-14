"""
Write out a shapefile for the grid edges.

Run this after changing the grid.

Output written to "derived" subdirectory.
"""

from stompy.model.delft import dfm_grid

g=dfm_grid.DFMGrid('sfei_v19_net.nc')

##

g.write_edges_shp('derived/grid-edges.shp')
