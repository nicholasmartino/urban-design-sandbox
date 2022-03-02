import os
import sys
from variables import *
import geopandas as gpd
import requests


def import_github(user='nicholasmartino', repository='morphology', branch='master', file='main.py'):
	r = requests.get(f"https://raw.githubusercontent.com/{user}/{repository}/{branch}/{file}")
	with open(f"models/{file}", "wb") as f:
		f.write(r.content)
	return

import_github(repository='morphology', file='ShapeTools.py')
import_github(repository='city', file='Network.py')

# Local settings
GRID_GDF = gpd.read_file(f'data/geojson/{GRID_FILE}')

# Server settings
if sys.platform == 'win32':
	SERVER_DIR = f'S:/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage'
else:
	SERVER_DIR = f'/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage'
# GRID_DIR = f'{SERVER_DIR}/Data/gmm_grids'
# GRID_FILES = ['dunbar_baseline.geojson', 'kitsilano_baseline.geojson']
# OUT_DIR = f'{SERVER_DIR}/Sandbox/shp/City-Wide'

if 'streets_utm_cov.feather' in os.listdir('data/feather'):
	STREETS = gpd.read_feather(f'data/feather/streets_utm_cov.feather')
else:
	STREETS = gpd.read_file(f'{SERVER_DIR}/Sandbox/shp/COV/Data/COV/streets_utm_cov.shp')
	STREETS.to_feather('data/feather/streets_utm_cov.feather')

TILE_GDF = gpd.read_feather('data/feather/all_tiles.feather')
NET_GDF = gpd.read_feather('data/feather/all_tiles_network.feather')
