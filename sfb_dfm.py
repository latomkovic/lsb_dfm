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

log=logging.getLogger('sfb_dfm')

import dredge_grid

DAY=np.timedelta64(86400,'s') # useful for adjusting times

## --------------------------------------------------

# Parameters to control more specific aspects of the run
if 0: # nice short setup for testing:
    run_name="short_20120801_p16" 
    run_start=np.datetime64('2012-08-01')
    run_stop=np.datetime64('2012-08-02')
if 1: # wy2013 with spinup
    run_name="wy2013" 
    run_start=np.datetime64('2012-08-01')
    run_stop=np.datetime64('2013-10-01')


# debugging - set all volumetric flow rates to 1m3/s.
ALL_FLOWS_UNIT=False

## --------------------------------------------------

# Derived parameters used in multiple places

# base_dir=os.path.dirname(__file__) # for command line invocation
base_dir="." # during development

run_base_dir=os.path.join(base_dir,'runs',run_name)

# path prefix for items in inputs-static
static_dir=os.path.relpath(os.path.join(base_dir,'inputs-static'),
                           run_base_dir)

# reference date - can only be specified to day precision, so 
# truncate to day precision (rounds down)
ref_date=run_start.astype('datetime64[D]')

net_file = os.path.join(base_dir,'sfei_v19_net.nc')

# No longer using any new-style boundary conditions
old_bc_fn = os.path.join(run_base_dir,'FlowFMold_bnd.ext')

obs_shp_fn = os.path.join(base_dir,'inputs-static','observation-points.shp')

dredge_depth=-0.5 # m NAVD88, depth to enforce at inflows and discharges

##

# Utility methods
def add_suffix_to_feature(feat,suffix):
    name=feat[0]
    suffize=lambda s: s.replace(name,name+suffix)
    feat_suffix=[suffize(feat[0]), feat[1]] # points stay the same
    if len(feat)==3: # includes names for nodes
        feat_suffix.append( [suffize(s) for s in feat[2]] )
    return feat_suffix



## 
# Make sure run directory exists:
os.path.exists(run_base_dir) or os.makedirs(run_base_dir)

# and that the bc files are clear:
for fn in [old_bc_fn]:
    os.path.exists(fn) and os.unlink(fn)

##

# Load the grid now -- it's used for clarifying some inputs, but
# is also modified to deepen areas near inflows, before being written
# out near the end of the script
grid=dfm_grid.DFMGrid(net_file)
    
## 

if 1:
    if 1:# fetch wind data, write to format supported by DFM
        target_filename_base=os.path.join(run_base_dir,"wind") # will have amu, amv appended

        wind_u_fn=target_filename_base+".amu"
        wind_v_fn=target_filename_base+".amv"
        if (os.path.exists(wind_u_fn) and
            os.path.exists(wind_v_fn)):
            log.info('Wind files already exist')
        else:
            data_start=run_start-1*DAY
            data_stop=run_stop+1*DAY

            url='http://sfbaynutrients.sfei.org/erddap/griddap/wind_ludwig_20170621'
            ds=xr.open_dataset(url)

            dio.dataset_to_dfm_wind(ds,data_start,data_stop,target_filename_base,
                                    extra_header="# downloaded from %s"%url)
    if 1: # and add wind to the boundary forcing
        wind_stanza=["QUANTITY=windx",
                     "FILENAME=wind.amu",
                     "FILETYPE=4",
                     "METHOD=2",
                     "OPERAND=O",
                     "",
                     "QUANTITY=windy",
                     "FILENAME=wind.amv",
                     "FILETYPE=4",
                     "METHOD=2",
                     "OPERAND=O",
                     "\n"]
        with open(old_bc_fn,'at') as fp:
            fp.write("\n".join(wind_stanza))

##

# features which have manually set locations for this grid
adjusted_pli_fn = os.path.join(base_dir,'nudged_features.pli')

if 1: #  Freshwater flows from git submodule
    adjusted_features=dio.read_pli(adjusted_pli_fn)
    # Add the freshwater flows - could come from erddap, but use github submodule
    # for better control on version

    # create a pair of bc and pli files, each including all the sources.
    # exact placement will
    # be done by hand in the GUI
    freshwater_dir=os.path.join(base_dir, 'sfbay_freshwater')

    full_flows_ds = xr.open_dataset(os.path.join(freshwater_dir, 'outputs', 'sfbay_freshwater.nc'))

    # period of the full dataset which will be include for this run
    sel=(full_flows_ds.time > run_start - 5*DAY) & (full_flows_ds.time < run_stop + 5*DAY)
    flows_ds = full_flows_ds.isel(time=sel)

    for stni in range(len(flows_ds.station)):
        stn_ds=flows_ds.isel(station=stni)

        src_name=stn_ds.station.item() # kind of a pain to get scalar values back out...

        # At least through the GUI, pli files must have more than one node.
        # Don't get too big for our britches, just stick a second node 50m east
        # if the incoming data is a point, but check for manually set locations
        # in adjusted_features
        if 1: #-- Write a PLI file
            feat=(src_name,
                  np.array( [[stn_ds.utm_x,stn_ds.utm_y],
                             [stn_ds.utm_x + 50.0,stn_ds.utm_y]] ))
            # Scan adjusted features for a match to use instead
            for adj_feat in adjusted_features:
                if adj_feat[0] == src_name:
                    feat=adj_feat
                    break
            # Write copies for flow, salinity and temperatures
            for suffix in ['_flow','_salt','_temp']:
                # function to add suffix
                feat_suffix=add_suffix_to_feature(feat,suffix)
                pli_fn=os.path.join(run_base_dir,"%s%s.pli"%(src_name,suffix))
                dio.write_pli(pli_fn,[feat_suffix])

            dredge_grid.dredge_boundary(grid,feat[1],dredge_depth)

        if 1: #-- Write the time series and stanza in FlowFM_bnd.ext
            df=stn_ds.to_dataframe().reset_index()
            df['elapsed_minutes']=(df.time.values - ref_date)/np.timedelta64(60,'s')
            df['salinity']=0*df.flow_cms
            df['temperature']=20+0*df.flow_cms

            if ALL_FLOWS_UNIT:
                df['flow_cms']=1.0+0*df.flow_cms
                
            for quantity,suffix in [ ('dischargebnd','_flow'),
                                     ('salinitybnd','_salt'),
                                     ('temperaturebnd','_temp') ]:
                lines=['QUANTITY=%s'%quantity,
                       'FILENAME=%s%s.pli'%(src_name,suffix),
                       'FILETYPE=9',
                       'METHOD=3',
                       'OPERAND=O',
                       ""]
                with open(old_bc_fn,'at') as fp:
                    fp.write("\n".join(lines))

                # read the pli back to know how to name the per-node timeseries
                feats=dio.read_pli(os.path.join(run_base_dir,
                                                "%s%s.pli"%(src_name,suffix)))
                feat=feats[0] # just one polyline in the file

                if len(feat)==3:
                    node_names=feat[2]
                else:
                    node_names=[""]*len(feat[1])
                    
                for node_idx,node_name in enumerate(node_names):
                    # if no node names are known, create the default name of <feature name>_0001
                    if not node_name:
                        node_name="%s%s_%04d"%(src_name,suffix,1+node_idx)

                    tim_fn=os.path.join(run_base_dir,node_name+".tim")
                    
                    columns=['elapsed_minutes']
                    if quantity=='dischargebnd':
                        columns.append('flow_cms')
                    elif quantity=='salinitybnd':
                        columns.append('salinity')
                    elif quantity=='temperaturebnd':
                        columns.append('temperature')
                    
                    df.to_csv(tim_fn, sep=' ', index=False, header=False, columns=columns)

##

# POTW inputs:
# The new-style boundary inputs file (FlowFM_bnd_new.ext) cannot represent
# sources and sinks, so these come in via the old-style file.

if 1:
    potws=xr.open_dataset(os.path.join(base_dir,'sfbay_potw','outputs','sfbay_delta_potw.nc'))
    adjusted_features=dio.read_pli(adjusted_pli_fn)

    # select a time subset of the flow data, starting just before the
    # simulation period, and running beyond the end:
    time_pnts = np.searchsorted(potws.time, [run_start-DAY,run_stop+DAY])
    time_pnts = time_pnts.clip(0,len(potws.time)-1)
    time_idxs=range(time_pnts[0],time_pnts[1]) # enumerate them for loops below

    with open(old_bc_fn,'at') as fp:
        for site in potws.site.values:
            # NB: site is bytes at this point
            potw=potws.sel(site=site)
            try:
                site_s=site.decode()
            except AttributeError:
                site_s=site

            if site_s in ['false_sac','false_sj']:
                print("(skip %s) "%site_s,end="")
                continue

            if potw.utm_x.values.mean() > 610200:
                # Delta POTWs are in this file, too, but not in this
                # grid.  Luckily they are easy to identify based on
                # x coordinate.
                print("(skip %s -- too far east) "%site_s,end="")
                continue
            
            print("%s "%site_s,end="")

            fp.write( ("QUANTITY=discharge_salinity_temperature_sorsin\n"
                       "FILENAME=%s.pli\n"
                       "FILETYPE=9\n"
                       "METHOD=1\n"
                       "OPERAND=O\n"
                       "AREA=0 # no momentum\n"
                       "\n")%site_s )

            # Write the location - writing a single point appears to work,
            # based on how it shows up in the GUI.  Otherwise we'd have to
            # manufacture a point outside the domain.
            with open(os.path.join(run_base_dir,'%s.pli'%site_s),'wt') as pli_fp:
                # Scan adjusted features for a match to use instead
                # This is handled slightly differently with POTWs - use the

                # put the depth at -50, should be at the bed
                feat=[site_s,
                      np.array([[potw.utm_x.values,potw.utm_y.values,-50.0]]),
                      ['']]

                for adj_feat in adjusted_features:
                    if adj_feat[0] == site_s:
                        # Merge them if the adjusted feature is more than 10 m away
                        # (to allow for some rounding in the ascii round-trip.)
                        offset=utils.dist( adj_feat[1][-1][:2] - feat[1][-1][:2] )
                        if offset > 10.0:
                            # Just add on the extra point - but may have to promote one 
                            # or the other to 3D.
                            old_geo=feat[1]
                            new_geo=adj_feat[1][-1:]
                            if old_geo.shape[1] != new_geo.shape[1]:
                                if old_geo.shape[1]<3:
                                    old_geo=np.concatenate( (old_geo,0*old_geo[:,:1]), axis=1)
                                else:
                                    # copy depth from old_geo
                                    new_geo=np.concatenate( (new_geo,
                                                             old_geo[-1,-1]+0*new_geo[:,:1]),
                                                            axis=1)

                            # if the original feature was outside the grid, then all is well,
                            # and it will show up in the GUI as a line from the original location
                            # outside the domain to the new location in the domain.
                            if grid.select_cells_nearest(old_geo[-1,:2],inside=True) is None:
                                feat[1]=np.concatenate( (old_geo,new_geo),axis=0 )
                                if len(feat)==3: # includes labels, but they don't matter here, right?
                                    feat[2].append('')
                            else:
                                # but if the original location is inside the grid, this will be interpreted
                                # as a sink-source pair, so we instead just put the single, adjusted
                                # location in.  This is done after potentially copying z-coordinate
                                # data from the original location.
                                feat[1]=new_geo
                        break

                dio.write_pli(pli_fp,[feat])

                dredge_grid.dredge_discharge(grid,feat[1],dredge_depth)

            with open(os.path.join(run_base_dir,'%s.tim'%site_s),'wt') as tim_fp:
                for tidx in time_idxs:
                    tstamp_minutes = (potw.time[tidx]-ref_date) / np.timedelta64(1,'m')

                    if ALL_FLOWS_UNIT:
                        flow_cms=1.0
                    else:
                        flow_cms=potw.flow[tidx]
                        
                    tim_fp.write("%g %g %g %g\n"%(tstamp_minutes,
                                                  flow_cms,
                                                  0.0, # salinity
                                                  20.0)) # temperature...


##

# Delta boundary conditions
# copied Silvia's pli files to inputs-static
# even though there are 14 of these, one per node of the sea boundary,
# they all appear to have the same data, even the tidal height time series.
# Sea_temp: probably grabbed from Point Reyes?
# Sea_sal: constant 33
# Sea_0001.pli - 

if 1: # Fetch river USGS river flows, add to FlowFM_bnd.ext:
    # Per Silvia's Thesis:
    # Jersey: Discharge boundary affected by tides, discharge and temperature taken
    # from USGS 11337190 SAN JOAQUIN R A JERSEY POINT, 0 salinity
    # (Note that Dutch Slough should probably be added in here)
    # Rio Vista: 11455420 SACRAMENTO A RIO VISTA, temperature from DWR station RIV.
    # 0 salinity.
    
    if 1: 
        # Cache the original data from USGS, then clean it and write to DFM format
        jersey_raw_fn=os.path.join(run_base_dir,'jersey-raw.nc')
        if not os.path.exists(jersey_raw_fn):
            jersey_raw=usgs_nwis.nwis_dataset(station="11337190",
                                              start_date=run_start,end_date=run_stop,
                                              products=[60, # "Discharge, cubic feet per second"
                                                        10], # "Temperature, water, degrees Celsius"
                                              days_per_request=30)
            jersey_raw.to_netcdf(jersey_raw_fn,engine='scipy')

        rio_vista_raw_fn=os.path.join(run_base_dir,'rio_vista-raw.nc')
        if not os.path.exists(rio_vista_raw_fn):
            rio_vista_raw=usgs_nwis.nwis_dataset(station="11455420",
                                                 start_date=run_start,end_date=run_stop,
                                                 products=[60, # "Discharge, cubic feet per second"
                                                           10], # "Temperature, water, degrees Celsius"
                                                 days_per_request=30)
            rio_vista_raw.to_netcdf(rio_vista_raw_fn,engine='scipy')


    if 1: # clean the data
        jersey_raw=xr.open_dataset(jersey_raw_fn)
        rio_vista_raw=xr.open_dataset(rio_vista_raw_fn)

        # fill any missing data via linear interpolation
        def fill_data(da):
            valid=np.isfinite(da.values)
            valid_times=da.time[valid]
            dt=np.median(np.diff(da.time))
            dt_gap = np.diff(da.time[valid]).max()
            if dt_gap > dt:
                log.warning("%s: gaps up to %.1f minutes"%(da.name, dt_gap/np.timedelta64(60,'s')))
            da.values = utils.fill_invalid(da.values)    

        fill_data(jersey_raw.stream_flow_mean_daily)
        fill_data(jersey_raw.temperature_water)

        fill_data(rio_vista_raw.stream_flow_mean_daily)
        fill_data(rio_vista_raw.temperature_water)

    if 1: # Write it all out
        for src_name,source in [ ('Jersey',jersey_raw),
                                 ('RioVista',rio_vista_raw)]:
            src_feat=dio.read_pli(os.path.join(base_dir,'inputs-static','%s.pli'%src_name))[0]
            
            df=source.to_dataframe().reset_index()
            df['elapsed_minutes']=(df.time.values - ref_date)/np.timedelta64(60,'s')
            
            if ALL_FLOWS_UNIT:
                df['flow_cms']=np.ones_like(df.elapsed_minutes)
            else:
                df['flow_cms']=0.028316847 * np.ones_like(df.elapsed_minutes)
                
            df['salinity']   = 0*np.ones_like(df.elapsed_minutes)
            df['temperature']=20*np.ones_like(df.elapsed_minutes)

            dredge_grid.dredge_boundary(grid,src_feat[1],dredge_depth)
            
            # Add stanzas to FlowFMold_bnd.ext:
            for quant,suffix in [('dischargebnd','_flow'),
                                 ('salinitybnd','_salt'),
                                 ('temperaturebnd','_temp')]:
                with open(old_bc_fn,'at') as fp:
                    lines=["QUANTITY=%s"%quant,
                           "FILENAME=%s%s.pli"%(src_name,suffix),
                           "FILETYPE=9",
                           "METHOD=3",
                           "OPERAND=O",
                           ""]
                    fp.write("\n".join(lines))
                    
                feat_suffix=add_suffix_to_feature(src_feat,suffix)
                dio.write_pli(os.path.join(run_base_dir,'%s%s.pli'%(src_name,suffix)),
                              [feat_suffix])

                # Write the data:
                if quant=='dischargebnd':
                    columns=['elapsed_minutes','flow_cms']
                elif quant=='salinitybnd':
                    columns=['elapsed_minutes','salinity']
                elif quant=='temperaturebnd':
                    columns=['elapsed_minutes','temperature']
                    
                if len(feat_suffix)==3:
                    node_names=feat_suffix[2]
                else:
                    node_names=[""]*len(feat_suffix[1])
                    
                for node_idx,node_name in enumerate(node_names):
                    # if no node names are known, create the default name of <feature name>_0001
                    if not node_name:
                        node_name="%s%s_%04d"%(src_name,suffix,1+node_idx)

                    tim_fn=os.path.join(run_base_dir,node_name+".tim")
                    df.to_csv(tim_fn, sep=' ', index=False, header=False, columns=columns)

## 

# Ocean:
#  Water level data from station 46214 (apparently from Yi Chao's ROMS?)
#    no spatial variation
#  Maybe salinity from Yi Chao ROMS?  That's what the thesis says, but the
#  actual inputs look like constant 33

if 1: # Ocean BCs from Point Reyes
    ptreyes_raw_fn=os.path.join(run_base_dir,'ptreyes-raw.nc')
    if 1:
        if not os.path.exists(ptreyes_raw_fn):
            ptreyes_raw=noaa_coops.coops_dataset("9415020",run_start,run_stop,
                                                 ["water_level","water_temperature"],
                                                 days_per_request=30)

            ptreyes_raw.to_netcdf(ptreyes_raw_fn,engine='scipy')

    if 1:
        # Clean that up, fabricate salinity
        ptreyes=xr.open_dataset(ptreyes_raw_fn).isel(station=0)

        # fill any missing data via linear interpolation
        def fill_data(da):
            valid=np.isfinite(da.values)
            valid_times=da.time[valid]
            dt=np.median(np.diff(da.time))
            dt_gap = np.diff(da.time[valid]).max()
            if dt_gap > dt:
                log.warning("%s: gaps up to %.1f minutes"%(da.name, dt_gap/np.timedelta64(60,'s')))
            da.values = utils.fill_invalid(da.values)    

        fill_data(ptreyes.water_level)
        fill_data(ptreyes.water_temperature)

        if ALL_FLOWS_UNIT:
            print("-=-=-=- USING 35 PPT WHILE TESTING! -=-=-=-")
            ptreyes['salinity']=35 + 0*ptreyes.water_temperature
        else:
            ptreyes['salinity']=33 + 0*ptreyes.water_temperature
            
    if 1: # Write it all out
        # Add a stanza to FlowFMold_bnd.ext:
        src_name='Sea'
        source=ptreyes
        df=source.to_dataframe().reset_index()
        df['elapsed_minutes']=(df.time.values - ref_date)/np.timedelta64(60,'s')
        
        forcing_data=[('waterlevelbnd','water_level','_ssh'),
                      ('salinitybnd','salinity','_salt'),
                      ('temperaturebnd','water_temperature','_temp')]

        src_feat=dio.read_pli(os.path.join(base_dir,'inputs-static','%s.pli'%src_name))[0]

        for quant,column,suffix in forcing_data:
            with open(old_bc_fn,'at') as fp:
                lines=["QUANTITY=%s"%quant,
                       "FILENAME=%s%s.pli"%(src_name,suffix),
                       "FILETYPE=9",
                       "METHOD=3",
                       "OPERAND=O",
                       ""]
                fp.write("\n".join(lines))

            feat_suffix=add_suffix_to_feature(src_feat,suffix)
            dio.write_pli(os.path.join(run_base_dir,'%s%s.pli'%(src_name,suffix)),
                          [feat_suffix])

            # Write the data:
            columns=['elapsed_minutes',column]

            if len(feat_suffix)==3:
                node_names=feat_suffix[2]
            else:
                node_names=[""]*len(feat_suffix[1])

            for node_idx,node_name in enumerate(node_names):
                # if no node names are known, create the default name of <feature name>_0001
                if not node_name:
                    node_name="%s%s_%04d"%(src_name,suffix,1+node_idx)

                tim_fn=os.path.join(run_base_dir,node_name+".tim")
                df.to_csv(tim_fn, sep=' ', index=False, header=False, columns=columns)

##

# Spatial salinity initial condition and friction
if 1:
    lines=[]
    if not ALL_FLOWS_UNIT: # real initial condition:
        lines+=[ "QUANTITY=initialsalinity",
                 "FILENAME=%s/saltopini.xyz"%static_dir,
                 "FILETYPE=7",
                 "METHOD=5",
                 "OPERAND=O",
                 ""]
    else: #  constant 35 ppt initial condition:
        print("-=-=-=- USING 35 PPT WHILE TESTING! -=-=-=-")
        lines+=[ "QUANTITY=initialsalinity",
                 "FILENAME=constant_35ppt.xyz",
                 "FILETYPE=7",
                 "METHOD=5",
                 "OPERAND=O",
                 ""]
        orig_salt=np.loadtxt('inputs-static/saltopini.xyz')
        orig_salt[:,2]=35
        np.savetxt(os.path.join(run_base_dir,'constant_35ppt.xyz'),
                   orig_salt,
                   delimiter=' ')
        
    lines+=["QUANTITY=frictioncoefficient",
            "FILENAME=%s/friction12e.xyz"%static_dir,
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

mdu['geometry','LandBoundaryFile'] = os.path.join(static_dir,"deltabay.ldb")

if 1:  # Copy grid file into run directory and update mdu
    mdu['geometry','NetFile'] = os.path.basename(net_file)
    dest=os.path.join(run_base_dir, net_file)
    if 0: # back when this script made no changes to the grid:
        shutil.copyfile(net_file, dest)
    else: # write out the modified grid
        dfm_grid.write_dfm(grid,dest,overwrite=True)

if 1: # fixed weir file is just referenced as static input
    mdu['geometry','FixedWeirFile'] = os.path.join(static_dir,'SBlevees_tdk.pli')

if 1: # set dates
    # RefDate can only be specified to day precision
    mdu['time','RefDate'] = utils.to_datetime(ref_date).strftime('%Y%m%d')
    mdu['time','Tunit'] = 'M' # minutes.  kind of weird, but stick with what was used already
    mdu['time','TStart'] = 0
    mdu['time','TStop'] = int( (run_stop - run_start) / np.timedelta64(1,'m') )

if 1: # update location of the boundary conditions
    # this has the source/sinks which cannot be written in the new style file
    mdu['external forcing','ExtForceFile']=os.path.basename(old_bc_fn)

if 0: #
    # Would be adding evaporation as negative rain here.
    pass

if 1: # output locations
    mdu['output','CrsFile'] = os.path.join(static_dir,"SB-observationcrosssection.pli")

##
if 1:
    # Observation points taken from shapefile for easier editing/comparisons in GIS
    # OLD mdu['output','ObsFile'] = os.path.join(static_dir,"Southbay-withoutDelta_obs.xyn")
    
    obs_pnts=wkb2shp.shp2geom(obs_shp_fn)
    obs_fn='observation_pnts.xyn'
    
    with open(os.path.join(run_base_dir,obs_fn),'wt') as fp:
        for idx,row in enumerate(obs_pnts):
            xy=np.array(row['geom'])
            fp.write("%12g %12g '%s'\n"%(xy[0], xy[1], row['name']))
    mdu['output','ObsFile'] = obs_fn

    if run_name.startswith('short'):
        mdu['output','MapInterval'] = 3600
    
## 
mdu.write(os.path.join(run_base_dir,run_name+".mdu"))

