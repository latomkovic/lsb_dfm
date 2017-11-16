import stompy.model.delft.io as dio
from stompy.spatial import wkb2shp
from shapely import geometry

## 

weirs=dio.read_pli('fixed_weirs-v02.pli')


geoms=[geometry.LineString(w[1][:,:2])
       for w in weirs]

##

wkb2shp.wkb2shp('fixed_weirs-v02.shp',geoms)
