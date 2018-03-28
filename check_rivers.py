# Compare sfbay_freshwater flows to gages.

import matplotlib.pyplot as  plt
import pandas as pd
import xarray as xr

from sfb_dfm_utils import common
from stompy import utils, filters

from stompy.io.local import usgs_nwis

## 

ds=xr.open_dataset('sfbay_freshwater/outputs/sfbay_freshwater.nc')

ds['dnum']=('time',),utils.to_dnum(ds.time)

## 

# Override BAHM with gage data when available
# what are possible overrides?
#  - COYOTE => 11172175 COYOTE C AB HWY 237 A MILPITAS CA
#  - SCLARAVCc => 11169025 GUADALUPE R ABV HWY 101 A SAN JOSE CA
#  - UALEMADAg => 11180700 ALAMEDA C FLOOD CHANNEL A UNION CITY CA
#  - USANLORZ => 11181040 SAN LORENZO C A SAN LORENZO CA
#  - MARINS3 => 11460000 CORTE MADERA C A ROSS CA
#  - MARINN => 11459500 NOVATO C A NOVATO CA
#  - NAPA => 11458000 NAPA R NR NAPA CA
# These are not flow stations:
#  - EBAY Cc4 => 374336122095801 SAN LEANDRO C A ALVARADO ST A SAN LEANDRO CA
#  - PETALUMA => 381519122385601 PETALUMA R NR PETALUMA CA

name_usgs_bahm= [ ('Coyote','11172175','COYOTE'),
                  ('Alviso','11169025','SCLARAVCc'),
                  ('Alameda Flood','11180700','UALAMEDA'),
                  ('San Lorenzo Ck','11181040', 'USANLORZ'), #  SAN LORENZO C A SAN LORENZO CA
                  ('Corte Madera','11460000','MARINS3'), # CORTE MADERA C A ROSS CA
                  ('Novato Ck','11459500','MARINN'), # NOVATO C A NOVATO CA
                  ('Napa River','11458000','NAPA'), # NAPA R NR NAPA CA
]

usgs_gages=[ usgs for name,usgs,bahm in name_usgs_bahm]


period_start=np.datetime64('2015-10-01')
period_end  =np.datetime64('2017-01-01')

gages=usgs_nwis.nwis_dataset_collection(usgs_gages,
                                        products=[60],
                                        start_date=period_start,
                                        end_date=period_end,
                                        days_per_request='M',
                                        cache_dir=common.cache_dir)

##

# Downsample to daily
gage_lp=gages.copy()
del gage_lp['tz_cd']
del gage_lp['datenum']

df=gage_lp.to_dataframe().unstack('site')
gage_daily=df.resample('D').mean()
# no need for second level index on columns
gage_daily=gage_daily.loc[:,('stream_flow_mean_daily')]
##

plt.figure(1).clf()
fig,axs=plt.subplots(len(name_usgs_bahm),1,num=1,sharex=True)

gage_dnum=utils.to_dnum(gage_daily.index.values)

for ax,(name,usgs,bahm) in zip(axs,name_usgs_bahm):
    #ax.plot(gages.datenum,gages.stream_flow_mean_daily.sel(site=usgs),
    #        marker=',',
    #        label='Gage: %s'%name)
    ax.plot(gage_dnum,
            gage_daily.loc[:,usgs],
            marker=',',
            label='Gage: %s'%name)

    ax.plot(ds.dnum, ds.flow_cfs.sel(station=bahm), label="BAHM: %s"%name)
    ax.legend()
    ax.axis(ymin=-5,ymax=gages.stream_flow_mean_daily.sel(site=usgs).max())

axs[0].axis(xmin=gage_dnum[0],xmax=gage_dnum[-1],
            ymin=-5,ymax=140)

axs[0].xaxis.axis_date()

# fig.tight_layout()
# fig.savefig('river_check_00.png')

## 

# For now, only worry about Alviso and Alameda Flood Control,
# wrap the steps into a concise method

def nudge_by_gage(ds,period_start,period_end,usgs_station,station,decorr_days):
    usgs_gage=usgs_nwis.nwis_dataset(usgs_station,
                                     products=[60],
                                     start_date=period_start,
                                     end_date=period_end,
                                     days_per_request='M',
                                     cache_dir=common.cache_dir)

    # Downsample to daily
    df=usgs_gage['stream_flow_mean_daily'].to_dataframe()
    df.index=df.index.levels[0] # comes in MultiIndex even though it's a single level
    df_daily=df.resample('D').mean()

    # Get the subset of BAHM data which overlaps this gage data
    time_slc=slice(np.searchsorted( ds.time, df_daily.index.values[0]),
                   1+np.searchsorted( ds.time, df_daily.index.values[-1] ) )

    bahm_subset=ds.sel(station=station).isel( time=time_slc )

    assert len(bahm_subset.time) == len(df_daily),"Maybe BAHM data doesn't cover entire period"

    errors=bahm_subset.flow_cfs - df_daily.stream_flow_mean_daily 

    # Easiest: interpolate errors over nans, apply to bahm data array.
    # the decorrelation time is tricky, though.

    # Specify a decorrelation time scale then relax from error to zero
    # over that period
    valid=np.isfinite(errors)
    errors_interp=np.interp( utils.to_dnum(ds.time),
                             utils.to_dnum(df_daily.index[valid]),
                             errors[valid])
    all_valid=np.zeros( len(ds.time),'f8' )
    all_valid[time_slc]=1*valid

    weights=(2*filters.lowpass_fir(all_valid,decorr_days)).clip(0,1)
    weighted_errors=weights*errors_interp

    # Does this work? 
    subset=dict(station=station)

    cfs_vals=ds.flow_cfs.loc[subset] - weighted_errors
    ds.flow_cfs.loc[subset] = cfs_vals.clip(0,np.inf) 
    ds.flow_cms.loc[subset] = 0.028316847 * ds.flow_cfs.loc[subset]

    # user feedback
    cfs_shifts=weighted_errors[time_slc]
    print("Shift in CFS: %.2f +- %.2f"%(np.mean(cfs_shifts),
                                        np.std(cfs_shifts)))

period_start=np.datetime64('2015-10-01')
period_end  =np.datetime64('2017-01-01')

nudge_by_gage(ds,period_start,period_end,'11169025',station='SCLARAVCc',decorr_days=20)
nudge_by_gage(ds,period_start,period_end,'11180700',station='UALAMEDA',decorr_days=20)

