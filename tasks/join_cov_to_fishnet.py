import geopandas as gpd
from variables import *


fishnet_file = f"{SERVER_DIR}/SSHRC Partnership Engage/Sandbox/shp/City-Wide/Fishnet_CoV_gmm_r4_Tile_Street_Types.shp"
cov_file = f"{SERVER_DIR}/02_Background info/CoV_received20220301_2050plans/Georeferenced/20220228 Neighbourhoods 2.shp"

fishnet = gpd.read_file(fishnet_file)
cov = gpd.read_file(cov_file)

for layer in cov['Layer'].unique():
    cov_layer = cov[cov['Layer'] == layer].copy()
    cov_layer[layer] = 1
    join = gpd.sjoin(fishnet, cov_layer.loc[:, [layer, 'geometry']])
    fishnet.loc[join.index, layer] = join[layer]
    fishnet[layer] = fishnet[layer].fillna(0)

cov_uu = cov.dissolve()
for i in range(5):
    cov_uu['geometry'] = cov_uu.buffer(-150)
    join2 = gpd.sjoin(fishnet, cov_uu.loc[:, ['geometry']])
    fishnet.loc[join2.index, 'level'] = i + 1

fishnet.to_file(
    f"{SERVER_DIR}/SSHRC Partnership Engage/Sandbox/shp/City-Wide/Fishnet_CoV_gmm_r4_Tile_Street_Types_CoV_Areas.shp"
)
