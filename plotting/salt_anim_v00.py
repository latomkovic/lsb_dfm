import glob
import os

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import collections
import numpy as np
import xarray as xr

from stompy.grid import unstructured_grid
from stompy.model.delft import dfm_grid
from stompy.plot import plot_wkb
from stompy import utils
import stompy.plot.cmap as scmap
from stompy.spatial import wkb2shp

##


##

if 0:
    # The long-ish run for Mark
    # But it looks like many of the levees aren't doing anything?
    # or maybe it's just that all control structures are wide open
    run_dir="../runs/medium_winter2016_01"
    out_dir="../runs/medium_winter2016_01/DFM_OUTPUT_medium_winter2016_01"
if 0:
    # Short, high-res winter run
    run_dir="../runs/short_winter2016_04"
    out_dir=os.path.join(run_dir,"DFM_OUTPUT_short_winter2016_04")
if 1:
    # as yet incomplete, high-res summer run
    run_dir="../runs/short_summer2016_00"
    out_dir=os.path.join(run_dir,"DFM_OUTPUT_short_summer2016_00")

if 0:
    # as yet incomplete, high-res summer run
    run_dir="../runs/short_winter2016_05"
    out_dir=os.path.join(run_dir,"DFM_OUTPUT_short_winter2016_05")

assert os.path.exists(out_dir)
## 

fig_dir=os.path.join(run_dir,'figures')

maps=glob.glob(os.path.join(out_dir,'*_map.nc'))

dss=[xr.open_dataset(m) for m in maps]

##

os.path.exists(fig_dir) or os.makedirs(fig_dir)


## 
G=dfm_grid.DFMGrid(os.path.join(run_dir,'lsb_v99_net.nc'))

grids=[ dfm_grid.DFMGrid(ds)
        for ds in dss ]

##

# Get some leveeish things to make it clearer what's going on.
lin_feats=wkb2shp.shp2geom("../fixed_weirs/linear_features-v02.shp")

# Trim down to ones that overlap our domain:
domain_poly=G.boundary_polygon()

from shapely.prepared import prep
# optimize for tests, and get rid of features that are part of the
# shoreline anyway.
domain_poly_prep=prep(domain_poly.buffer(-15))

# down to 5k features
sel_lin_feats=[ lf
                for lf in lin_feats
                if ( domain_poly_prep.intersects(lf['geom'])
                     and lf['model_name']=='')]

# merge to a single set of segments:
lin_segs=[ np.array(lf['geom'])
           for lf in sel_lin_feats ]


##

# Figure out the mapping of subdomain to global:
global_elems=np.zeros(G.Ncells(),'int32')-1

for proc,g in enumerate(grids):
    print(proc)
    centers=g.cells_centroid()
    for c in range(g.Ncells()):
        Gcell=G.select_cells_nearest( centers[c],inside=True )
        assert Gcell>=0
        # Gcell is the 0-based cell index in G match with local
        # 0-based cell index c.
        # local 0-based cell index is in turn known by this
        # global element nr:
        global_elems[Gcell]=dss[proc].FlowElemGlobalNr.values[c]
G.add_cell_field('global_elem',global_elems)
# and establish the mapping from global element to cell index in G
# fill with invalid first
global_to_G=np.zeros(G.cells['global_elem'].max()+1,'int32') -1
# then update valid
global_to_G[global_elems]=np.arange(G.Ncells())

# Can we just reorder the original grid cells?
# reorder accepts an argsort-like index array.
# the first element of order gives the original index which should now be 1st.
# I think this is what is in global_elems?
# I want whatever element has FlowElemGlobalNr==0 to be specified in order.
# okay - so it's not global_elems, it's the permutation of that -
# so global_to_G[0] = 56280.
# ah - but global_elems is 1 based?

# What I want is to wrap the collection of map outputs into a single object
# I can treat like a merged map output.
# so I would do whatever selection of variables

##

def merge(selector):
    Gvar=None

    for ds in dss:
        per_proc=selector(ds)
        if Gvar is None:
            Gvar=np.zeros(G.Ncells(),per_proc.dtype)
        mapping=global_to_G[ds.FlowElemGlobalNr.values]
        Gvar[mapping] = per_proc
    return Gvar

##
zoom=(575471.08620748692, 595825.67607978475, 4138312.4467484825, 4153604.563194477)
zoom_in=(576406.34727196465, 595359.95406865189, 4138583.1679984042, 4152869.1084242631)

def setup_figure(num=1):
    plt.figure(num).clf()
    fig,ax=plt.subplots(num=num)
    fig.set_size_inches([6.4,4.8],forward=True)
    ax.set_position([0,0,1,1])
    ax.xaxis.set_visible(0)
    ax.yaxis.set_visible(0)
    ax.set_facecolor('0.7')

    lcoll=collections.LineCollection(lin_segs,lw=0.8,color='k')
    ax.add_collection(lcoll)
    plot_wkb.plot_wkb(domain_poly,facecolor='none',lw=0.8,ax=ax)
    plot_wkb.plot_wkb(domain_poly,edgecolor='none',facecolor='w',lw=0.8,ax=ax,zorder=-.5)

    cax=fig.add_axes([0.08,0.2,0.3,0.03])

    return fig,ax,cax

##

fig,ax,cax=setup_figure(1)

Gvar=np.zeros(G.Ncells(),'f8')

coll=G.plot_cells(values=Gvar,clip=zoom,ax=ax)
cmap=scmap.load_gradient('StepSeq25.cpt')
plt.setp(coll,lw=1,edgecolor='face',clim=[0,33],cmap=cmap)

date_txt=ax.text(0.08,0.27,"1111-22-33 44:55",
                 transform=ax.transAxes)

cell_mask=G.cell_clip_mask(zoom)
def set_time(t_idx):
    Gvar[:] = merge(lambda ds: ds['sa1'].isel(time=t_idx).isel(laydim=0).values)
    depth = merge(lambda ds: ds['waterdepth'].isel(time=t_idx).values)
    coll.set_array( np.ma.array(Gvar[cell_mask],mask=depth[cell_mask]<0.05 ) )
    t=dss[0].time.values[t_idx]
    date_txt.set_text( utils.to_datetime(t).strftime('%Y-%m-%d %H:%M') )

plt.colorbar(coll,cax=cax,orientation='horizontal',label='Salinity (ppt)')
ax.axis(zoom_in)

set_time(0) 

##

frame_dir=os.path.join(fig_dir,'salt_frames_v00')
os.path.exists(frame_dir) or os.makedirs(frame_dir)

for frame,t_idx in enumerate(range(0,len(dss[0].time))):
    set_time(t_idx)
    #fig.canvas.draw()
    #plt.pause(0.05)
    img_fn=os.path.join(frame_dir,'salt_%04d.png'%frame)
    print(img_fn)
    fig.savefig(img_fn,dpi=100)

##
if 0:
    # Close up on A3W
    zoom_a3w=(585208.24208633008, 586676.51929676568, 4143107.0491332593, 4144208.257041086)
    ax.axis(zoom_a3w)
    coll.set_clim([14,20])
    set_time(199)
    #set_time(196)
    plt.draw()
    fig.savefig("guad-a3w-salt-map.png",dpi=100)

    ##

    # Stratification
    fig,ax,cax=setup_figure(2)

    Gvar=np.zeros(G.Ncells(),'f8')

    coll=G.plot_cells(values=Gvar,clip=zoom,ax=ax)
    plt.setp(coll,lw=1,edgecolor='face',clim=[-1,1],cmap='seismic')
    plt.colorbar(coll,cax=cax,orientation='horizontal',label='Stratification (ppt/m)',
                 ticks=[-1,-0.5,0,0.5,1])

    date_txt=ax.text(0.08,0.27,"1111-22-33 44:55",
                     transform=ax.transAxes)

    cell_mask=G.cell_clip_mask(zoom)
    def set_time(t_idx): 
        s_surf=merge(lambda ds: ds.sa1.isel(time=t_idx).isel(laydim=0).values)
        s_bed=merge(lambda ds: ds.sa1.isel(time=t_idx).isel(laydim=-1).values)
        depth = merge(lambda ds: ds['waterdepth'].isel(time=t_idx).values)

        Gvar[:] = (s_surf-s_bed)/depth.clip(1,np.inf)

        coll.set_array( np.ma.array(Gvar[cell_mask],mask=depth[cell_mask]<0.05 ) )
        t=dss[0].time.values[t_idx]
        date_txt.set_text( utils.to_datetime(t).strftime('%Y-%m-%d %H:%M') )

    ax.axis(zoom_in)

    set_time(5) 

    ##

    frame_dir=os.path.join(fig_dir,'strat_frames_v00')
    os.path.exists(frame_dir) or os.makedirs(frame_dir)

    for frame,t_idx in enumerate(range(0,len(dss[0].time))):
        set_time(t_idx)
        img_fn=os.path.join(frame_dir,'strat_%04d.png'%frame)
        print(img_fn)
        fig.savefig(img_fn,dpi=100)
