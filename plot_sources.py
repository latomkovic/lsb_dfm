"""
Debugging the sources -- early checks do not make it clear whether the sources are
properly included or not.
"""
import os

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from stompy.grid import unstructured_grid
from stompy.model.delft import dfm_grid

## 

output_dir="runs/short_20120801_p24/DFM_OUTPUT_short_20120801_p24"

Nproc=24
procs=range(Nproc)

all_ds=[xr.open_dataset(os.path.join(output_dir,'short_20120801_p24_%04d_20120801_000000_map.nc'%proc))
        for proc in procs]

all_g=[unstructured_grid.UnstructuredGrid.from_ugrid(ds)
       for ds in all_ds]
##
plt.figure(1).clf()
fig,ax=plt.subplots(num=1)


# Plot lowest layer salinity for all grids:
all_coll=[]
for proc in procs:
    salt0=all_ds[proc].sa1.isel(time=0,laydim=9).values
    saltN=all_ds[proc].sa1.isel(time=-1,laydim=9).values
    dsalt=saltN-salt0
    coll=all_g[proc].plot_cells(ax=ax,lw=0.5,values=dsalt)
    coll.set_clim([-.1,0])
    all_coll.append(coll)

plt.setp(all_coll,edgecolor='face')
plt.colorbar(all_coll[0])
