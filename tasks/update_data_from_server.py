import geopandas as gpd
import sys


# Server settings
if sys.platform == 'win32':
	SERVER_DIR = f'S:/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage'
else:
	SERVER_DIR = f'/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage'
DIR = f'{SERVER_DIR}/Sandbox/shp/elementslab/Version_3'
TILE_GDF = gpd.read_file(f'{DIR}/all_tiles.shp').to_crs(26910)
# TREE_GDF = gpd.read_file('/Volumes/Macintosh HD/Users/nicholasmartino/Google Drive/elementslab/sandbox_tiles_trees.geojson')
NET_GDF = gpd.read_file(f'{DIR}/Tiles_Network.shp').to_crs(26910)
TILES = []

TILE_GDF.to_feather('../data/feather/all_tiles.feather')
NET_GDF.to_feather('../data/feather/all_tiles_network.feather')
