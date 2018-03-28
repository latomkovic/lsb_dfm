# Compare sfbay_freshwater flows to gages.

import matplotlib.pyplot as  plt
import pandas as pd
import xarray as xr

from sfb_dfm_utils import common
from stompy import utils

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
# This is not a flow station:
#  - EBAY Cc4 => 374336122095801 SAN LEANDRO C A ALVARADO ST A SAN LEANDRO CA
#  - MARINS3 => 11460000 CORTE MADERA C A ROSS CA
#  - MARINN => 11459500 NOVATO C A NOVATO CA
# Not a flow station:
#  - PETALUMA => 381519122385601 PETALUMA R NR PETALUMA CA
#  - NAPA => 11458000 NAPA R NR NAPA CA

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

daily=pd.PeriodIndex(start=period_start,end=period_end,freq='D')
quarters=pd.PeriodIndex(start=period_start,end=period_end,freq='15M')


## 

# Downsample to daily
gage_lp=gages.copy()
del gage_lp['tz_cd']

#for site in gage_lp.site.values:
#    gage_lp.stream_flow_mean_daily.values[:] = 

##

plt.figure(1).clf()
fig,axs=plt.subplots(len(name_usgs_bahm),1,num=1,sharex=True)

for ax,(name,usgs,bahm) in zip(axs,name_usgs_bahm):
    ax.plot(gages.datenum,gages.stream_flow_mean_daily.sel(site=usgs),
            marker=',',
            label='Gage: %s'%name)

    ax.plot(ds.dnum, ds.flow_cfs.sel(station=bahm), label="BAHM: %s"%name)
    ax.legend()
    ax.axis(ymin=-5,ymax=gages.stream_flow_mean_daily.sel(site=usgs).max())

axs[0].axis(xmin=gages.datenum[0],xmax=gages.datenum[-1],ymin=-5,ymax=140)

axs[0].xaxis.axis_date()

# fig.tight_layout()
# fig.savefig('river_check_00.png')

## 



