# Exploring how to improve the initial salinity field in
# LSB runs.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from stompy.grid import unstructured_grid
from stompy.model.delft import dfm_grid
from stompy.spatial import wkb2shp
import stompy.plot.cmap as scmap
from stompy import utils

## 

run_start=np.datetime64('2016-06-01')

saltopini_xyz=np.loadtxt('inputs-static/saltopini.xyz')

##

g=dfm_grid.DFMGrid("lsb_v99_bathy_net.nc")

##

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

g.plot_edges(ax=ax,lw=0.5,color='k')

scat=ax.scatter(saltopini_xyz[:,0],
                saltopini_xyz[:,1],
                40,
                saltopini_xyz[:,2],
                cmap=scmap.load_gradient('hot_desaturated.cpt'))
plt.colorbar(scat)

##

# Get some observations:

# First, from a recent USGS cruise:
from stompy.io.local import usgs_sfbay

usgs_data_end=np.datetime64('2016-04-28')

usgs_pad=np.timedelta64(30,'D')

usgs_target=run_start

# lame, but if we have to, reach to a prior year
while usgs_target + usgs_pad > usgs_data_end:
    usgs_target -= np.timedelta64(365,'D')

usgs_cruises=usgs_sfbay.cruise_dataset(usgs_target - usgs_pad,
                                       usgs_target + usgs_pad )

##

# lame filling
salt3d=usgs_cruises['salinity']
salt2d=salt3d.mean(dim='prof_sample')
assert salt2d.dims[0]=='date'
salt2d_fill=utils.fill_invalid(salt2d.values,axis=0)

from scipy.interpolate import interp1d

salt_f=interp1d(utils.to_dnum(salt2d.date.values),
                salt2d_fill,
                axis=0,bounds_error=False)(utils.to_dnum(usgs_target))

usgs_init_salt=np.c_[salt2d.x.values,salt2d.y.values,salt_f]
##


plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

g.plot_edges(ax=ax,lw=0.5,color='k',zorder=-1)

scat=ax.scatter(usgs_init_salt[:,0],
                usgs_init_salt[:,1],
                40,
                usgs_init_salt[:,2],
                cmap=scmap.load_gradient('hot_desaturated.cpt'))
plt.colorbar(scat)

##

# And pull some SFEI data:

mooring_xy=[]
mooring_salt=[]

L2_dir='/opt/data/sfei/moored_sensors_csv/L2/'

# tuples (<name in observation points shapefile>, <L2 data file name> )
sfei_moorings=[
    ('ALV',"ALV_all_data_L2.csv"),
    ('SFEI_Coyote',"COY_all_data_L2.csv"),
    ('DB',"DMB_all_data_L2.csv"),
    ('SFEI_Guadalupe',"GL_all_data_L2.csv"),
    ('SFEI_Mowry',"MOW_all_data_L2.csv"),
    ('SFEI_Newark',"NW_all_data_L2.csv"),
    ('SFEI_A8Notch',"POND_all_data_L2.csv"),
    ('SMB',"SM_all_data_L2.csv")
]

# lat/lon from observation-points
obs_shp=wkb2shp.shp2geom("inputs-static/observation-points.shp")


for name,l2_file in sfei_moorings:
    print(name)
    sfei=pd.read_csv(os.path.join(L2_dir,l2_file),
                     parse_dates=['Datetime'])
    sfei_salt=sfei['S_PSU']
    valid=~(sfei_salt.isnull())
    sfei_salt_now=np.interp(utils.to_dnum(run_start),
                            utils.to_dnum(sfei.Datetime[valid]),sfei_salt[valid])
    geom=obs_shp['geom'][ np.nonzero(obs_shp['name']==name)[0][0] ]
    xy=np.array(geom)
    if np.isfinite(sfei_salt_now):
        mooring_xy.append(xy)
        mooring_salt.append(sfei_salt_now)

##

xy=np.array(mooring_xy)

sfei_init_salt=np.c_[xy[:,0],xy[:,1],mooring_salt]

##
        
init_salt=np.concatenate( (usgs_init_salt,
                           sfei_init_salt) )

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

g.plot_edges(ax=ax,lw=0.5,color='k',zorder=-1)

scat=ax.scatter(init_salt[:,0],
                init_salt[:,1],
                40,
                init_salt[:,2],
                cmap=scmap.load_gradient('hot_desaturated.cpt'))

plt.colorbar(scat)

##

from stompy.model import unstructured_diffuser

differ=unstructured_diffuser.Diffuser(grid=g)

##

for x,y,salt in init_salt:
    try:
        differ.set_dirichlet(salt,xy=[x,y],on_duplicate='skip')
    except differ.PointOutsideDomain as exc:
        print(exc)
        continue
        

##
differ.construct_linear_system()
differ.solve_linear_system(animate=False)

##

coll2=g.plot_cells(values=differ.C_solved,
                   cmap=scmap.load_gradient('hot_desaturated.cpt'))
coll2.set_clim( scat.get_clim() )

##

cc=g.cells_centroid()

##

cc_salt=np.concatenate( ( cc, differ.C_solved[:,None] ),axis=1 )

##

# Because DFM is going to use some interpolation, and will not reach outside
# the convex hull, we have to be extra cautious and throw some points out farther
# afield.

xys_orig=np.loadtxt('inputs-static/orig-saltopini.xyz')

combined_xys=np.concatenate( (cc_salt,xys_orig), axis=0 )

## 
np.savetxt('inputs-static/lsb_saltopini.xyz',combined_xys)
