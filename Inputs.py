import geopandas as gpd
import sys

# Local settings
GRID_DIR = 'data/grids'
GRID_FILES = [
	'dunbar_baseline.geojson',
    'kitsilano_baseline.geojson',
	# 'main_baseline.geojson',
	# 'marpole_baseline.geojson',
	# 'sunset_baseline.geojson'
]
OUT_DIR = 'data/sandboxes'

# Server settings
if sys.platform == 'win32':
	SERVER_DIR = f'S:/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage'
else:
	SERVER_DIR = f'/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage'
GRID_DIR = f'{SERVER_DIR}/Data/gmm_grids'
GRID_FILES = ['dunbar_baseline.geojson', 'kitsilano_baseline.geojson']
OUT_DIR = f'{SERVER_DIR}/Sandbox/shp/City-Wide'

STREETS = gpd.read_file(f'{SERVER_DIR}/Sandbox/shp/COV/Data/COV/streets_utm_cov.shp')
