# Compare sfbay_freshwater flows to gages.

import matplotlib.pyplot as  plt
import xarray as xr

from sfb_dfm_utils import common
from stompy import utils

from stompy.io.local import usgs_nwis

## 

ds=xr.open_dataset('sfbay_freshwater/outputs/sfbay_freshwater.nc')

ds['dnum']=('time',),utils.to_dnum(ds.time)

## 

coyote='11172175' # 
alviso='11169025' # GUADALUPE R ABV HWY 101
alameda='11180700'


gages=usgs_nwis.nwis_dataset_collection([coyote,alviso,alameda],
                                        products=[60],
                                        start_date=np.datetime64('2015-10-01'),
                                        end_date=np.datetime64('2017-01-01'),
                                        days_per_request='M',
                                        cache_dir=common.cache_dir)



## 

name_usgs_bahm= [ ('Coyote','11172175','COYOTE'),
                  ('Alviso','11169025','SCLARAVCc'),
                  ('Alameda Flood','11180700','UALAMEDA') ]
         
gages['dnum']=('time',),utils.to_dnum(gages.time.values)

plt.figure(1).clf()
fig,axs=plt.subplots(len(name_usgs_bahm),1,num=1,sharex=True)

for ax,(name,usgs,bahm) in zip(axs,name_usgs_bahm):
    ax.plot(gages.dnum,gages.stream_flow_mean_daily.sel(site=usgs),
            label='Gage: %s'%name)

    ax.plot(ds.dnum, ds.flow_cfs.sel(station=bahm), label="BAHM: %s"%name)
    ax.legend()
    ax.axis(ymin=-5,ymax=gages.stream_flow_mean_daily.sel(site=usgs).max())

axs[0].axis(xmin=gages.dnum[0],xmax=gages.dnum[-1],ymin=-5,ymax=140)

axs[0].xaxis.axis_date()

fig.tight_layout()

fig.savefig('river_check_00.png')

## 



