#!/usr/bin/env python
"""
Script to configure and generate mdu ready for running dflow fm.

Uses sfei_v19 grid, a Southbay-enhancement of the SF Bay/Delta community
model grid, with some areas deepened (LSB bathy), trimmed (Coyote Creek)
and dredge (see dredge_grid.py in this directory)
"""
import os
import pdb
import io
import shutil
import numpy as np
import logging
import xarray as xr
import six

from stompy import utils
import stompy.model.delft.io as dio
from stompy.model.delft import dfm_grid
from stompy.spatial import wkb2shp
from stompy.io.local import usgs_nwis,noaa_coops

if __name__=='__main__':
    logging.basicConfig(level=logging.INFO)

log=logging.getLogger('lsb_dfm')

import sfb_dfm_utils 
import set_bathy

DAY=np.timedelta64(86400,'s') # useful for adjusting times

## --------------------------------------------------

# Parameters to control more specific aspects of the run
if 1: # nice short setup for testing:
    run_name="short_10d_p16" 
    run_start=np.datetime64('2016-06-01')
    run_stop=run_start+10*DAY
if 0: # wy2013 with spinup
    run_name="wy2013" 
    run_start=np.datetime64('2012-08-01')
    run_stop=np.datetime64('2013-10-01')


ALL_FLOWS_UNIT=False # for debug, set all volumetric flow rates to 1m3/s if True

## --------------------------------------------------

# Derived parameters used in multiple places

# base_dir=os.path.dirname(__file__) # for command line invocation
base_dir="." # during development

run_base_dir=os.path.join(base_dir,'runs',run_name)

# real location of static directory
abs_static_dir=os.path.join(base_dir,'inputs-static')
# and how to refer to it relative to the run directory
rel_static_dir=os.path.relpath(abs_static_dir,
                               run_base_dir)


# reference date - can only be specified to day precision, so 
# truncate to day precision (rounds down)
ref_date=run_start.astype('datetime64[D]')

net_file = os.path.join(base_dir,'lsb_v99_net.nc')

# No longer using any new-style boundary conditions
old_bc_fn = os.path.join(run_base_dir,'FlowFMold_bnd.ext')

obs_shp_fn = os.path.join(abs_static_dir, 'observation-points.shp')

dredge_depth=-0.5 # m NAVD88, depth to enforce at inflows and discharges

##

# Make sure run directory exists:
os.path.exists(run_base_dir) or os.makedirs(run_base_dir)

# clear any stale bc files:
for fn in [old_bc_fn]:
    os.path.exists(fn) and os.unlink(fn)

##

# Load the grid now -- it's used for clarifying some inputs, but
# is also modified to deepen areas near inflows, before being written
# out near the end of the script
# Bathy for LSB is evolving, so a few steps required to possibly update it here.

net_bathy_file=net_file.replace("_net.nc","_bathy_net.nc")
assert net_bathy_file!=net_file
if ( (not os.path.exists(net_bathy_file))
     or (os.stat(net_file).st_mtime >= os.stat(net_bathy_file).st_mtime) ):
    grid=dfm_grid.DFMGrid(net_file)
    log.info("Will update bathymetry in grid - 2 minutes?")
    set_bathy.set_lsb_bathy(grid)
    log.info("Writing updated grid/bathy")
    dfm_grid.write_dfm(grid,net_bathy_file,overwrite=True)

# Either way, read that back in 
log.info("Reading grid with bathy")
grid=dfm_grid.DFMGrid(net_bathy_file)

## 

# Write out a shapefile for the grid edges.
# Run this after changing the grid.
# Output written to "derived" subdirectory.
edges_shp="derived/grid-edges.shp"

if ( (not os.path.exists(edges_shp))
     or (os.stat(net_file).st_mtime >= os.stat(edges_shp).st_mtime) ):
    log.info("Writing new edge shapefile")
    grid.write_edges_shp('derived/grid-edges.shp')
else:
    log.info("Edge shapefile exists.  Will leave as is.")

##

# WIND
sfb_dfm_utils.add_erddap_ludwig_wind(run_base_dir,
                                     run_start,run_stop,
                                     old_bc_fn)
                                         
##

# features which have manually set locations for this grid
adjusted_pli_fn = os.path.join(base_dir,'nudged_features.pli')
# pretty slow...
sfb_dfm_utils.add_sfbay_freshwater(run_base_dir,
                                   run_start,run_stop,ref_date,
                                   adjusted_pli_fn,
                                   freshwater_dir=os.path.join(base_dir, 'sfbay_freshwater'),
                                   grid=grid,
                                   dredge_depth=dredge_depth,
                                   old_bc_fn=old_bc_fn,
                                   all_flows_unit=ALL_FLOWS_UNIT)
                     
##

# POTW inputs:
# The new-style boundary inputs file (FlowFM_bnd_new.ext) cannot represent
# sources and sinks, so these come in via the old-style file.
potw_dir=os.path.join(base_dir,'sfbay_potw')
# why so slow?
sfb_dfm_utils.add_sfbay_potw(run_base_dir,
                             run_start,run_stop,ref_date,
                             potw_dir,
                             adjusted_pli_fn,
                             grid,dredge_depth,
                             old_bc_fn,
                             all_flows_unit=ALL_FLOWS_UNIT)

##

# Delta boundary conditions
# For short period, this has condition number issues
sfb_dfm_utils.add_delta_inflow(run_base_dir,
                               run_start,run_stop,ref_date,
                               static_dir=abs_static_dir,
                               grid=grid,dredge_depth=dredge_depth,
                               old_bc_fn=old_bc_fn,
                               all_flows_unit=ALL_FLOWS_UNIT)
##

# For the short run, no data for temperature
sfb_dfm_utils.add_ocean(run_base_dir,
                        run_start,run_stop,ref_date,
                        static_dir=abs_static_dir,
                        grid=grid,
                        old_bc_fn=old_bc_fn,
                        all_flows_unit=ALL_FLOWS_UNIT)

##

sfb_dfm_utils.add_initial_salinity(run_base_dir,
                                   static_dir=abs_static_dir,
                                   old_bc_fn=old_bc_fn,
                                   all_flows_unit=ALL_FLOWS_UNIT)
## 
if 1:            
    lines=["QUANTITY=frictioncoefficient",
           "FILENAME=%s/friction12e.xyz"%rel_static_dir,
           "FILETYPE=7",
           "METHOD=5",
           "OPERAND=O",
           ""]
    with open(old_bc_fn,'at') as fp:
        fp.write("\n".join(lines))

## --------------------------------------------------------------------------------
# Edits to the template mdu:
# 

mdu=dio.MDUFile('template.mdu')

mdu['geometry','LandBoundaryFile'] = os.path.join(rel_static_dir,"deltabay.ldb")

# Start with 2D:
# mdu['geometry','Kmx']=1
# On to 3D, for better or worse
mdu['geometry','Kmx']=10

if 1:  # Copy grid file into run directory and update mdu
    mdu['geometry','NetFile'] = os.path.basename(net_file)
    dest=os.path.join(run_base_dir, net_file)
    # write out the modified grid
    dfm_grid.write_dfm(grid,dest,overwrite=True)

if 1: # fixed weir file is just referenced as static input
    mdu['geometry','FixedWeirFile'] = os.path.join(rel_static_dir,'FlowFM_fxw.pli')

if 1: # set dates
    # RefDate can only be specified to day precision
    mdu['time','RefDate'] = utils.to_datetime(ref_date).strftime('%Y%m%d')
    mdu['time','Tunit'] = 'M' # minutes.  kind of weird, but stick with what was used already
    mdu['time','TStart'] = 0
    mdu['time','TStop'] = int( (run_stop - run_start) / np.timedelta64(1,'m') )

if 1: # update location of the boundary conditions
    # this has the source/sinks which cannot be written in the new style file
    mdu['external forcing','ExtForceFile']=os.path.basename(old_bc_fn)

if 0: # Would be adding evaporation as negative rain here.
    pass

if 1: # output locations
    mdu['output','CrsFile'] = os.path.join(rel_static_dir,"SB-observationcrosssection.pli")

##
if 1:
    # Observation points taken from shapefile for easier editing/comparisons in GIS
    obs_pnts=wkb2shp.shp2geom(obs_shp_fn)
    obs_fn='observation_pnts.xyn'
    
    with open(os.path.join(run_base_dir,obs_fn),'wt') as fp:
        for idx,row in enumerate(obs_pnts):
            xy=np.array(row['geom'])
            fp.write("%12g %12g '%s'\n"%(xy[0], xy[1], row['name']))
    mdu['output','ObsFile'] = obs_fn

    if run_name.startswith('short'):
        mdu['output','MapInterval'] = 1800
        mdu['output','HisInterval'] = 900

## 
mdu.write(os.path.join(run_base_dir,run_name+".mdu"))

