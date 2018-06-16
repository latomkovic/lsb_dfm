#!/usr/bin/env python
"""
Script to configure and generate mdu ready for running dflow fm.

Uses sfei_v19 grid, a Southbay-enhancement of the SF Bay/Delta community
model grid, with some areas deepened (LSB bathy), trimmed (Coyote Creek)
and dredge (see dredge_grid.py in this directory)
"""
import os
import glob
import pdb
import io
import shutil
import subprocess
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

if 0: # slightly longer lsb summer run
    run_name="short_summer2016_04"
    run_start=np.datetime64('2016-06-01')
    run_stop=run_start+90*DAY

if 1: # re-run a short winter run
    run_name="short_winter2015_05"
    run_start=np.datetime64("2015-12-15")
    run_stop=np.datetime64("2016-01-30")

nprocs=4


# despite issues with DWAQ output, this remains a good version.
dfm_bin_dir="/home/rusty/src/dfm/r53925-opt/bin"

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

## --------------------------------------------------------------------------------
# Edits to the template mdu:
#

mdu=dio.MDUFile('template.mdu')

if 1: # set dates
    # RefDate can only be specified to day precision
    mdu['time','RefDate'] = utils.to_datetime(ref_date).strftime('%Y%m%d')
    mdu['time','Tunit'] = 'M' # minutes.  kind of weird, but stick with what was used already
    mdu['time','TStart'] = 0
    mdu['time','TStop'] = int( (run_stop - run_start) / np.timedelta64(1,'m') )

mdu['geometry','LandBoundaryFile'] = os.path.join(rel_static_dir,"deltabay.ldb")

mdu['geometry','Kmx']=10 # 10 layers

# update location of the boundary conditions
# this has the source/sinks which cannot be written in the new style file
mdu['external forcing','ExtForceFile']=os.path.basename(old_bc_fn)

# Load the grid now -- it's used for clarifying some inputs, but
# is also modified to deepen areas near inflows, before being written
# out near the end of the script
# Bathy for LSB is evolving, so a few steps required to possibly update it here.

net_bathy_file=net_file.replace("_net.nc","_bathy_net.nc")
bathy_file="inputs-static/merged_2m.tif"
assert net_bathy_file!=net_file

if ( (not os.path.exists(net_bathy_file))
     or (os.stat(net_file).st_mtime >= os.stat(net_bathy_file).st_mtime)
     or (os.stat(bathy_file).st_mtime >= os.stat(net_bathy_file).st_mtime)  ):
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

# features which have manually set locations for this grid
adjusted_pli_fn = os.path.join(base_dir,'nudged_features.pli')
# pretty slow...
sfb_dfm_utils.add_sfbay_freshwater(run_base_dir,
                                   run_start,run_stop,ref_date,
                                   adjusted_pli_fn,
                                   freshwater_dir=os.path.join(base_dir, 'sfbay_freshwater'),
                                   grid=grid,
                                   dredge_depth=dredge_depth,
                                   old_bc_fn=old_bc_fn)

##

# POTW inputs:
# The new-style boundary inputs file (FlowFM_bnd_new.ext) cannot represent
# sources and sinks, so these come in via the old-style file.
potw_dir=os.path.join(base_dir,'sfbay_potw')
# why so slow?
write_salt=mdu['physics','salinity']=='1'
write_temp=mdu['physics','temperature']=='1'

sfb_dfm_utils.add_sfbay_potw(run_base_dir,
                             run_start,run_stop,ref_date,
                             potw_dir,
                             adjusted_pli_fn,
                             grid,dredge_depth,
                             old_bc_fn,
                             write_temp=write_temp,
                             write_salt=write_salt)

##

# Delta boundary conditions
# For short period, this has condition number issues
sfb_dfm_utils.add_delta_inflow(run_base_dir,
                               run_start,run_stop,ref_date,
                               static_dir=abs_static_dir,
                               grid=grid,dredge_depth=dredge_depth,
                               old_bc_fn=old_bc_fn)
##

# For the short run, no data for temperature
sfb_dfm_utils.add_ocean(run_base_dir,
                        run_start,run_stop,ref_date,
                        static_dir=abs_static_dir,
                        grid=grid,
                        factor=0.901,
                        old_bc_fn=old_bc_fn)

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

if 1:  # Copy grid file into run directory and update mdu
    mdu['geometry','NetFile'] = os.path.basename(net_file)
    dest=os.path.join(run_base_dir, mdu['geometry','NetFile'])
    # write out the modified grid
    dfm_grid.write_dfm(grid,dest,overwrite=True)


## Initial salinity - do this once the grid has been written out
sfb_dfm_utils.add_initial_salinity_dyn(run_base_dir,
                                       abs_static_dir,
                                       mdu,
                                       run_start)

# WIND
ludwig_ok=sfb_dfm_utils.add_erddap_ludwig_wind(run_base_dir,
                                               run_start,run_stop,
                                               old_bc_fn)
if not ludwig_ok:
    const_ok=sfb_dfm_utils.add_constant_wind(run_base_dir,mdu,[4,0],run_start,run_stop)
    assert const_ok

##



fixed_weir_out=os.path.join(base_dir,'fixed_weirs','out')
if 1: # fixed weir file is just referenced as static input
    # mdu['geometry','FixedWeirFile'] = os.path.join(rel_static_dir,'FlowFM_fxw.pli')
    # updated with some features inside A5/7/8
    # since the code for this is now part of this repo, this isn't really static, so
    # copy it in.
    shutil.copyfile( os.path.join(fixed_weir_out,'fixed_weirs-v02.pli'),
                     os.path.join(run_base_dir,'fixed_weirs-v02.pli') )
    mdu['geometry','FixedWeirFile'] = 'fixed_weirs-v02.pli'

if 1: # add in gates, also derived in fixed_weirs
    # gate-specific inputs to the old-style inp file
    mdu['geometry','StructureFile'] = 'gates-v04.ini'
    shutil.copyfile( os.path.join(fixed_weir_out,'gates-v04.ini'),
                     os.path.join(run_base_dir,'gates-v04.ini') )
    for f in glob.glob( os.path.join(fixed_weir_out,'gate-*.pli') ):
        shutil.copyfile( f, os.path.join(run_base_dir, os.path.basename(f) ) )


if 0:  # this is creating a lot of high-salt spikes - omit for now
    sfb_dfm_utils.add_cimis_evap_precip(run_base_dir,mdu,scale_evap=1.0)

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

    if run_name.startswith('medium'):
        mdu['output','MapInterval'] = 3600
        mdu['output','HisInterval'] = 900

##
mdu_fn=os.path.join(run_base_dir,run_name+".mdu")
mdu.write(mdu_fn)

##

# As of r52184, explicitly built with metis support, partitioning can be done automatically
# from here.

cmd="%s/mpiexec -n %d %s/dflowfm --partition:ndomains=%d %s"%(dfm_bin_dir,nprocs,dfm_bin_dir,nprocs,
                                                              mdu['geometry','NetFile'])
pwd=os.getcwd()
try:
    os.chdir(run_base_dir)
    res=subprocess.call(cmd,shell=True)
finally:
    os.chdir(pwd)


# similar, but for the mdu:
cmd="%s/generate_parallel_mdu.sh %s %d 6"%(dfm_bin_dir,os.path.basename(mdu_fn),nprocs)
try:
    os.chdir(run_base_dir)
    res=subprocess.call(cmd,shell=True)
finally:
    os.chdir(pwd)


# 10 days at 0.5h => 42G
# 75 days at 1h => 150G

cmd="%s/mpiexec -n %d %s/dflowfm --autostartstop %s"%(dfm_bin_dir,nprocs,dfm_bin_dir,
                                                      os.path.basename(mdu_fn))
pwd=os.getcwd()
try:
    os.chdir(run_base_dir)
    res=subprocess.call(cmd,shell=True)
finally:
    os.chdir(pwd)
