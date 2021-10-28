import geopandas as gpd

# Local settings
GRID_DIR = 'data/grids'
GRID_FILES = ['dunbar_baseline.geojson', 'kitsilano_baseline.geojson', 'main_baseline.geojson',
              'marpole_baseline.geojson', 'sunset_baseline.geojson']
OUT_DIR = 'data/sandboxes'

# Server settings
SERVER_DIR = f'/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage'
GRID_DIR = f'{SERVER_DIR}/Data/gmm_grids'
GRID_FILES = ['fishnet_CoV_gmm_r4_built.shp']
OUT_DIR = f'{SERVER_DIR}/Sandbox/shp/City-Wide'

STREETS = gpd.read_file(
	'/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Sandbox/shp/COV/Data/COV/streets_utm_cov.shp')
