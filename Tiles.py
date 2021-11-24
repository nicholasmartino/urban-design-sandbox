import geopandas as gpd
from Tile import Tile
from Inputs import SERVER_DIR

DIR = f'{SERVER_DIR}/Sandbox/shp/elementslab/Version_3'
TILE_GDF = gpd.read_file(f'{DIR}/all_tiles.shp').to_crs(26910)
# TREE_GDF = gpd.read_file('/Volumes/Macintosh HD/Users/nicholasmartino/Google Drive/elementslab/sandbox_tiles_trees.geojson')
NET_GDF = gpd.read_file(f'{DIR}/Tiles_Network.shp').to_crs(26910)
TILES = []

assert 'cell_type' in TILE_GDF.columns
assert 'Subtype' in TILE_GDF.columns
assert 'bldgs' in TILE_GDF['Type'].unique()
assert 'prcls' in TILE_GDF['Type'].unique()
assert 'block' in TILE_GDF['Type'].unique()

for tp in TILE_GDF['cell_type'].unique():
	for st in TILE_GDF[TILE_GDF['cell_type'] == tp]['Subtype'].dropna().unique():
		tile = Tile(
			name=tp,
			buildings=TILE_GDF[(TILE_GDF['Type'] == 'bldgs') & (TILE_GDF['cell_type'] == tp) & (TILE_GDF['Subtype'] == st)],
			parcels=TILE_GDF[(TILE_GDF['Type'] == 'prcls') & (TILE_GDF['cell_type'] == tp) & (TILE_GDF['Subtype'] == st)],
			network=NET_GDF[(NET_GDF['cell_type'] == tp) & (NET_GDF['Subtype'] == st)],
			trees=TILE_GDF[(TILE_GDF['Type'] == 'trees') & (TILE_GDF['cell_type'] == tp) & (TILE_GDF['Subtype'] == st)],
			block=TILE_GDF[(TILE_GDF['Type'] == 'block') & (TILE_GDF['cell_type'] == tp) & (TILE_GDF['Subtype'] == st)],
			bound=TILE_GDF[(TILE_GDF['Type'] == 'bound') & (TILE_GDF['cell_type'] == tp) & (TILE_GDF['Subtype'] == st)],
			subtype=st)
		TILES.append(tile)
		if st in ['Standard', 'T-Shape', 'L-Middle', 'L-Edge']:
			TILES.append(Tile(tp, tile.buildings, tile.parcels, tile.network, tile.trees, tile.block, tile.bound, f'{st} Down').flip_vertical())
			TILES.append(Tile(tp, tile.buildings, tile.parcels, tile.network, tile.trees, tile.block, tile.bound, f'{st} Left').flip_horizontal())
			TILES.append(Tile(tp, tile.buildings, tile.parcels, tile.network, tile.trees, tile.block, tile.bound, f'{st} Down Left').flip_horizontal().flip_vertical())
