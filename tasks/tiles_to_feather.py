import geopandas as gpd
from variables import *

# Backup current files
current_network = gpd.read_feather(f"{FEATHER_DIR}/all_tiles_network.feather")
current_network.to_feather(f"{FEATHER_DIR}/all_tiles_network_backup.feather")

current_tiles = gpd.read_feather(f'{FEATHER_DIR}/all_tiles.feather')
current_tiles.to_feather(f'{FEATHER_DIR}/all_tiles_backup.feather')

# New files
tile_gdf = gpd.read_file(f"{TILE_DIR}/2022-03-01-NewTiles.shp")
tile_gdf.to_feather(f"{FEATHER_DIR}/all_tiles.feather")

network_gdf = gpd.read_file(f"{NETWORK_DIR}/Tiles_Network.shp")
network_gdf.to_feather(f"{FEATHER_DIR}/all_tiles_network.feather")
