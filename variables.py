GRID_DIR = 'data/grids'
GRID_FILE = 'broadway_plan_baseline.geojson'  # Default grid file
GRID_FILES = [
    'dunbar_baseline.geojson',
    'kitsilano_baseline.geojson',
    'main_baseline.geojson',
    'marpole_baseline.geojson',
    'sunset_baseline.geojson'
]
OUT_DIR = 'data/sandboxes'
TYPES = {
    2: 'Open_Low_Density',
    3: 'Mid_High_Street',
    4: 'Coarse_Grain',
    5: 'Moderate_Density',
    8: 'Treed_Large_Home',
    9: 'Dense_Nodal',
    10: 'Typical_Van_SF',
    11: 'Typical_Van_West_SF',
    12: 'Green_Open_Space',
    13: 'Dense_TOD',
    14: 'Dense_Industrial',
    15: 'Moderate_Industrial'
}
COLORS = {
    2: '8CD1E5',
    3: 'F4A259',
    4: 'CCCCCC',
    5: 'CE6A85',
    8: 'F4C26F',
    9: 'BC4B51',
    10: 'F4E285',
    11: 'CFD05A',
    12: '8CB369',
    13: '702D30',
    14: 'B37D80',
    15: '985277'
}
SERVER_DIR = "R:/Research/eLabs/50_projects/20_City_o_Vancouver"
TILE_DIR = f"{SERVER_DIR}/22_SSHRC Partnership Engage/02_shapefiles"
NETWORK_DIR = f"{SERVER_DIR}/SSHRC Partnership Engage/Sandbox/shp/elementslab/Version_3"
FEATHER_DIR = "../data/feather"
