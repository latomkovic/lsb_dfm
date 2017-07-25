#!/usr/bin/env python
"""
Script to configure and generate mdu ready for running dflow fm.

Uses sfei_v18 grid, a Southbay-enhancement of the SF Bay/Delta community
model grid.
"""
import os
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

log=logging.getLogger('sfb_dfm')
log.setLevel(logging.INFO)

DAY=np.timedelta64(86400,'s') # useful for adjusting times

## --------------------------------------------------

# Parameters to control more specific aspects of the run
if 0: # nice short setup for testing:
    run_name="short_20120801_p24" 
    run_start=np.datetime64('2012-08-01')
    run_stop=np.datetime64('2012-08-10')
if 1: # wy2013 with spinup
    run_name="wy2013" 
    run_start=np.datetime64('2012-08-01')
    run_stop=np.datetime64('2013-10-01')
    
## --------------------------------------------------

# Derived parameters used in multiple places

# base_dir=os.path.dirname(__file__) # for command line invocation
base_dir="." # during development

run_base_dir=os.path.join(base_dir,'runs',run_name)

# path prefix for items in inputs-static
static_dir=os.path.relpath(os.path.join(base_dir,'inputs-static'),
                           run_base_dir)

# reference date - can only be specified to day precision, so it may not
# truncate to day precision (rounds down)
ref_date=run_start.astype('datetime64[D]')

net_file = os.path.join(base_dir,'sfei_v18_net.nc')

new_bc_fn = os.path.join(run_base_dir,'FlowFMnew_bnd.ext')
old_bc_fn = os.path.join(run_base_dir,'FlowFMold_bnd.ext')

obs_shp_fn = os.path.join(base_dir,'inputs-static','observation-points.shp')

##

# Make sure run directory exists:
os.path.exists(run_base_dir) or os.makedirs(run_base_dir)

# and that the bc files are clear:
for fn in [new_bc_fn,old_bc_fn]:
    os.path.exists(fn) and os.unlink(fn)

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

if 0: # doesn't work for Delta Shell to write all the features into a single file
    fresh_pli_fn = os.path.join(run_base_dir, 'freshwater.pli')
else:
    fresh_pli_fn = None # will write features one-by-one

# It does work to have a single bc file holding the data and variables.
fresh_bc_fn = os.path.join(run_base_dir, 'freshwater.bc')

if 1: #  Freshwater flows from git submodule
    adjusted_features=dio.read_pli(adjusted_pli_fn)
    # Add the freshwater flows - could come from erddap, but use github submodule
    # for better control on version

    # create a pair of bc and pli files, each including all the sources.
    # exact placement will
    # be done by hand in the GUI
    freshwater_dir=os.path.join(base_dir, 'sfbay_freshwater')

    full_flows_ds = xr.open_dataset(os.path.join(freshwater_dir, 'outputs', 'sfbay_freshwater.nc'))

    sel=(full_flows_ds.time > run_start - 5*DAY) & (full_flows_ds.time < run_stop + 5*DAY)
    flows_ds = full_flows_ds.isel(time=sel)

    pli_data=[]

    bc_fp = open(fresh_bc_fn, 'wt')

    for stni in range(len(flows_ds.station)):
        stn_ds=flows_ds.isel(station=stni)

        src_name=stn_ds.station.item() # kind of a pain to get scalar values back out...

        # At least through the GUI, pli files must have more than one node.
        # Don't get too big for our britches, just stick a second node 50m east
        # if the incoming data is a point
        if 1: #-- Write a PLI file
            feat=(src_name,
                  np.array( [[stn_ds.utm_x,stn_ds.utm_y],
                             [stn_ds.utm_x + 50.0,stn_ds.utm_y]] ))
            # Scan adjusted features for a match to use instead
            for adj_feat in adjusted_features:
                if adj_feat[0] == src_name:
                    feat=adj_feat
                    break

            pli_data.append(feat)

        if 1: #-- Write a BC file
            df=stn_ds.to_dataframe().reset_index()

            df['unix_time']=utils.to_unix(df.time.values)

            # I think the 0001 has to be there, as it is used to
            # specify different values at different nodes of the pli
            # seems better to assume that incoming data is a daily average,
            # and keep it constant across the day
            # block-from means that the timestamp of a datum marks the beginning
            # of a period, and the value is held constant until the next timestamp
            # how about unix epoch for time units?
            bc_fp.write("[forcing]\n")
            # This Name needs to match the name in the pli
            bc_fp.write("Name               = %s_0001\n"%src_name)
            bc_fp.write("Function           = timeseries\n")
            bc_fp.write("Time-interpolation = block-from\n") # or linear, block-to 
            bc_fp.write("Quantity           = time\n")
            bc_fp.write("Unit               = seconds since 1970-01-01 00:00:00\n")
            bc_fp.write("Quantity           = dischargebnd\n")
            bc_fp.write("Unit               = m3 s-1\n")  # was "m³/s\n", but prefer simpler characters

            df.to_csv(bc_fp, sep=' ', index=False, header=False, columns=['unix_time', 'flow_cms'])

    bc_fp.close()
    if fresh_pli_fn is not None: # will write them all to a single file
        dio.write_pli(fresh_pli_fn,pli_data)
    else:
        for feat in pli_data:
            dio.write_pli( os.path.join(run_base_dir,feat[0]+".pli"),
                           [feat] )

    # Use the "new" format for FlowFM_bnd.ext, writing one file, with a stanza
    # for each boundary condition, all referencing the same forcingfile, but
    # with the features split out to individual pli files.
    with open(new_bc_fn,'at') as fp:
        for feat in pli_data:
            # This is how things were written out for hyperpycnal:
            # the same forcingfile can be specified multiple times, and it looks like
            # they will pull out the stanza matching
            # unclear from the manual whether multiple features can be included in a particular
            # .pli file here.
            # while the format
            pli_fn=feat[0]+".pli"
            fp.write("[boundary]\n")
            fp.write("quantity=dischargebnd\n")
            fp.write("locationfile=%s\n"%pli_fn)
            fp.write("forcingfile=%s\n"%os.path.basename(fresh_bc_fn))
            fp.write("return_time=0 # disable thatcher-harlemann\n")
            fp.write("\n")

##

# POTW inputs:
# The new-style boundary inputs file (FlowFM_bnd_new.ext) cannot represent
# sources and sinks, so these come in via a separate, old-style file.

if 1:
    grid=dfm_grid.DFMGrid(net_file)
    
    potws=xr.open_dataset(os.path.join(base_dir,'sfbay_potw','outputs','sfbay_potw.nc'))
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
            if six.PY3:
                site_s=site.decode()
            else:
                site_s=site

            if site_s in ['false_sac','false_sj']:
                print("(skip %s) "%site_s,end="")
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

            with open(os.path.join(run_base_dir,'%s.tim'%site_s),'wt') as tim_fp:
                for tidx in time_idxs:
                    tstamp_minutes = (potw.time[tidx]-ref_date) / np.timedelta64(1,'m')
                    tim_fp.write("%g %g %g %g\n"%(tstamp_minutes,
                                                  potw.flow[tidx],
                                                  0.0, # salinity
                                                  20.0)) # temperature...


##

# Ocean and Delta boundary conditions
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
            # Add a stanze to FlowFMnew_bnd.ext:
            with open(new_bc_fn,'at') as fp:
                for quant in ['dischargebnd','temperaturebnd']:
                    fp.write("\n".join(["\n[boundary]",
                                        "quantity=%s"%quant,
                                        "locationfile=%s.pli"%src_name,
                                        "forcingfile=%s.bc"%src_name,
                                        "return_time=0 # disable thatcher-harlemann",
                                        ""]))

            # Write the data:
            df=source.to_dataframe().reset_index()
            df['unix_time']=utils.to_unix(df.time.values)
            df['flow_cms']=0.028316847 * df['stream_flow_mean_daily']
            with open(os.path.join(run_base_dir,'%s.bc'%src_name),'wt') as bc_fp:
                for quant,units,column in [ ('dischargebnd','m3 s-1','flow_cms'),
                                            ('temperaturebnd','degC','temperature_water') ]:
                    bc_fp.write("\n".join( [ "[forcing]",
                                             # 'Name' needs to match the name in the pli
                                             "Name               = %s_0001"%src_name,
                                             "Function           = timeseries",
                                             "Time-interpolation = linear",
                                             "Quantity           = time",
                                             "Unit               = seconds since 1970-01-01 00:00:00",
                                             "Quantity           = %s"%quant,
                                             "Unit               = %s"%units,  # was "m³/s\n", but prefer simpler characters
                                             "\n"]))

                    df.to_csv(bc_fp, sep=' ', index=False, header=False, columns=['unix_time', column])

            # And copy the pli:
            jersey_pli   =dio.read_pli(os.path.join(base_dir,'inputs-static','Jersey.pli'))
            rio_vista_pli=dio.read_pli(os.path.join(base_dir,'inputs-static','RioVista.pli'))

            dio.write_pli(os.path.join(run_base_dir,'Jersey.pli'  ), jersey_pli)
            dio.write_pli(os.path.join(run_base_dir,'RioVista.pli'), rio_vista_pli)

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

        ptreyes['salinity']=33 + 0*ptreyes.water_temperature

    if 1: # Write it all out
        # Add a stanze to FlowFMnew_bnd.ext:
        src_name='Sea'
        source=ptreyes
        forcing_data=[('waterlevelbnd','m','water_level'),
                      ('salinitybnd','ppt','salinity'),
                      ('temperaturebnd','degC','water_temperature')]
        
        with open(new_bc_fn,'at') as fp:
            for quant,_,_ in forcing_data:
                fp.write("\n".join(["\n[boundary]",
                                    "quantity=%s"%quant,
                                    "locationfile=%s.pli"%src_name,
                                    "forcingfile=%s.bc"%src_name,
                                    "return_time=0 # disable thatcher-harlemann",
                                    ""]))

        # Write the data:
        df=source.to_dataframe().reset_index()
        df['unix_time']=utils.to_unix(df.time.values)

        with open(os.path.join(run_base_dir,'%s.bc'%src_name),'wt') as bc_fp:
            for quant,units,column in forcing_data:
                txt="\n".join(["[forcing]",
                               # 'Name' needs to match the name in the pli
                               "Name               = %s_0001"%src_name,
                               "Function           = timeseries",
                               "Time-interpolation = linear",
                               "Quantity           = time",
                               "Unit               = seconds since 1970-01-01 00:00:00",
                               "Quantity           = %s"%quant,
                               "Unit               = %s"%units,  # was "m³/s\n", but prefer simpler characters
                               ""])
                bc_fp.write(txt)
                df.to_csv(bc_fp, sep=' ', index=False, header=False,
                          columns=['unix_time', column])
                bc_fp.write("\n")

        # And copy the pli:
        sea_pli   =dio.read_pli(os.path.join(base_dir,'inputs-static','Sea.pli'))

        dio.write_pli(os.path.join(run_base_dir,'Sea.pli'), sea_pli)

##

# Spatial salinity initial condition and friction
if 1:
    txt="\n".join([ "QUANTITY=initialsalinity",
                    "FILENAME=%s/saltopini.xyz"%static_dir,
                    "FILETYPE=7",
                    "METHOD=5",
                    "OPERAND=O",
                    "",
                    "QUANTITY=frictioncoefficient",
                    "FILENAME=%s/friction12e.xyz"%static_dir,
                    "FILETYPE=7",
                    "METHOD=5",
                    "OPERAND=O",
                    ""])
    with open(old_bc_fn,'at') as fp:
        fp.write(txt)

## --------------------------------------------------------------------------------
# Edits to the template mdu:
# 

mdu=dio.MDUFile('template.mdu')

mdu['geometry','LandBoundaryFile'] = os.path.join(static_dir,"deltabay.ldb")

if 1:  # Copy grid file into run directory and update mdu
    shutil.copyfile(net_file, os.path.join(run_base_dir, net_file))
    mdu['geometry','NetFile'] = os.path.basename(net_file)

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
    # this has boundary flows
    mdu['external forcing','ExtForceFileNew']=os.path.basename(new_bc_fn)

if 0: #
    # Would be adding evaporation as negative rain here.
    pass



if 1: # output locations
    mdu['output','CrsFile'] = os.path.join(static_dir,"SB-observationcrosssection.pli")

    # Observation points taken from shapefile for easier editing/comparisons in GIS
    
    # mdu['output','ObsFile'] = os.path.join(static_dir,"Southbay-withoutDelta_obs.xyn")
##
if 1:    
    obs_pnts=wkb2shp.shp2geom(obs_shp_fn)
    obs_fn='observation_pnts.xyn'
    
    with open(os.path.join(run_base_dir,obs_fn),'wt') as fp:
        for idx,row in enumerate(obs_pnts):
            xy=np.array(row['geom'])
            fp.write("%12g %12g '%s'\n"%(xy[0], xy[1], row['name']))
    mdu['output','ObsFile'] = obs_fn
    
## 
mdu.write(os.path.join(run_base_dir,run_name+".mdu"))

