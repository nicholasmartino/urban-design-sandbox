import geopandas as gpd
from shapely.geometry import Polygon

gdf = gpd.read_file('C:/Users/nichmar.stu/Desktop/Georeferenced/cov_shapes.geojson')

for layer in gdf['Layer']:
    layer_gdf = gdf[gdf['Layer'] == layer]
    layer_gdf['geometry'] = [Polygon(geom) for geom in layer_gdf['geometry']]
    layer_gdf.to_file(f'C:/Users/nichmar.stu/Desktop/Georeferenced/CoV - {layer}.geojson', driver='GeoJSON')
