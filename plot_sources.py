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
    # steps 3 and 4 should be clear of the 3h ramp up.
    salt0=all_ds[proc].sa1.isel(time=0,laydim=9).values
    saltN=all_ds[proc].sa1.isel(time=4,laydim=9).values
    dsalt=saltN-salt0
    coll=all_g[proc].plot_cells(ax=ax,lw=0.5,values=dsalt)
    coll.set_clim([-.1,0])
    all_coll.append(coll)
    # show domain boundaries:
    e2c=all_g[proc].edge_to_cells()
    all_g[proc].plot_edges(ax=ax,lw=1.0,color='k',mask=e2c.min(axis=1)<0)

plt.setp(all_coll,edgecolor='face')
plt.colorbar(all_coll[0])

## 
import glob
for pli_fn in glob.glob(os.path.join(run_base_dir,'*.pli')):
    if '_salt' in pli_fn or '_temp' in pli_fn: continue
    feat=dio.read_pli(pli_fn)[0]
    name=feat[0].replace('_flow','')
    ax.text(feat[1][0,0],feat[1][0,1],name)
    ax.plot(feat[1][:,0],feat[1][:,1],'m-',lw=2)

##

# select a region based on the current axes:
zoom=ax.axis()

# iterate through grids:
time_idx=3 # 3600s per step
totalQ=0
for proc in procs:
    cells=np.nonzero( all_g[proc].cell_clip_mask(zoom) )[0]
    if len(cells)==0:
        continue

    times=[time_idx,time_idx+1]
    salt=all_ds[proc].sa1.isel(time=times,nFlowElem=cells).mean(dim='laydim').values
    fresh=(35-salt)/35.
    areas=all_g[proc].cells_area()[None,cells]
    # netcdf does not like doubling up on indexing
    depths=all_ds[proc].waterdepth.isel(time=times).values[:,cells]

    fresh_vol=(fresh*areas*depths).sum(axis=1)
    dfresh_vol=fresh_vol[1] - fresh_vol[0]
    
    print("freshwater volume per h, proc %d: %.3f"%(proc,dfresh_vol))
    totalQ+=dfresh_vol
print("freshwater volume per h, total: %.3f"%totalQ)


# Richardson Bay source is just 614.
# EBDA is spot on!
# Mud Slough is way low
# Coyote creek is way low.
# San Jose is spot on.

# Seems like there is a problem with the unit forcing for freshwater flows
# time series for Coyote is 1.0 for both nodes.
##

plt.figure(2).clf()
fig,ax=plt.subplots(num=2)

# Plot ssh for all grids:
all_coll=[]
for proc in procs:
    sshN=all_ds[proc].s1.isel(time=-1).values
    coll=all_g[proc].plot_cells(ax=ax,lw=0.5,values=sshN)
    # coll.set_clim([-.1,0])
    all_coll.append(coll)

plt.setp(all_coll,edgecolor='face',clim=[-2,2])
plt.colorbar(all_coll[0])
##

plt.figure(3).clf()
fig,ax=plt.subplots(num=3)

# Plot waterdepth for all grids:
all_coll=[]
for proc in procs:
    depth=all_ds[proc].waterdepth.isel(time=0).values
    coll=all_g[proc].plot_cells(ax=ax,lw=0.5,values=depth,cmap='jet')
    all_coll.append(coll)

plt.setp(all_coll,edgecolor='face',clim=[-0.1,10.0])
plt.colorbar(all_coll[0])

##

plt.figure(4).clf()
fig,ax=plt.subplots(num=4)

# proc 22, which has Alviso Slough
proc=22
cell_depth=all_ds[proc].FlowElem_bl.values
node_depth=all_ds[proc].NetNode_z.values

ccoll=all_g[proc].plot_cells(ax=ax,lw=0.5,values=cell_depth,cmap='jet',zorder=0)
ncoll=all_g[proc].plot_nodes(ax=ax,values=node_depth,cmap='jet',zorder=2,lw=1)
#ncoll.set_markeredgecolor('k')
plt.setp([ccoll,ncoll],edgecolor='face',clim=[-2,2]) 
from stompy.plot import plot_utils 
plot_utils.cbar(ccoll,extras=[ncoll])


##

# Is this a problem with dry cells not accepting flow?
pli=dio.read_pli("runs/short_20120801_p24/RioVista_flow.pli")

pnts=pli[0][1]
ax.plot(pnts[:,0],pnts[:,1],'g-o')

##

# Something related to which links are "boundary" links?

# nNetLink is 5096, also number of edges in grid.
# nFlowLink is 4240
# BndLink ranges from 1 to 5095
is_bnd=np.zeros(all_g[21].Nedges(),'b1')
is_bnd[ all_ds[21].BndLink.values -1 ]=True

all_g[21].plot_edges(mask=is_bnd,color='k',lw=2)

##

# now have flow data output in q1 ~ (time, nFlowLink, laydim)
# for proc 21:
# nNetLink: 5096
# nBndLink: 868 -- non-internal edges of the subdomain
# nFlowLink: 4240

# nNetLink - nFlowLink = 856
#  suggesting that 868 - 856 = 12 FlowLinks are on the boundary.

proc=22
ax.plot( all_ds[proc].FlowLink_xu.values,
         all_ds[proc].FlowLink_yu.values,
         'mo')
all_g[proc].plot_edges(color='k',lw=0.5)

## 
# Narrow that down to the velocity points on the boundary:
flow_xy=np.c_[all_ds[proc].FlowLink_xu.values,
              all_ds[proc].FlowLink_yu.values]
sac_flowlinks=np.nonzero(utils.within_2d(flow_xy, ax.axis()))[0]
# [4237, 4238, 4239]
sac_inflow=all_ds[proc].q1.isel(nFlowLink=[4237,4238,4239])
sac_flows = sac_inflow.sum(dim='laydim').sum(dim='nFlowLink')
# sac_inflow.isel(time=1).sum() => 0.51840
# or over an hour, that's 1866.2411988189333
