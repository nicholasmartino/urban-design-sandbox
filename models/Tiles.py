import gc

from models.Tile import Tile
from store import *

assert 'cell_type' in TILE_GDF.columns
assert 'Subtype' in TILE_GDF.columns
assert 'bldgs' in TILE_GDF['Type'].unique()
assert 'prcls' in TILE_GDF['Type'].unique()
assert 'block' in TILE_GDF['Type'].unique()

TILES = []
for tp in TILE_GDF['cell_type'].unique():
	unique = TILE_GDF[TILE_GDF['cell_type'] == tp]['Subtype'].unique()
	print(f"{tp}: {unique}")
	for st in unique:
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
		gc.collect()
