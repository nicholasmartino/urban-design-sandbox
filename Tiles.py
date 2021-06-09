import geopandas as gpd
from Tile import Tile

DIR = '/Volumes/ELabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Sandbox/shp/elementslab/Version_2'
TILE_GDF = gpd.read_file('data/sandbox_tiles.geojson')
TREE_GDF = gpd.read_file('/Volumes/Macintosh HD/Users/nicholasmartino/Google Drive/elementslab/sandbox_tiles_trees.geojson')
NET_GDF = gpd.read_file(f'{DIR}/Tiles_Network.shp')
# VARIATIONS = ['Standard', 'Vertical', 'L-Middle', 'L-Edge', 'T-Shape', 'Z-Shape']
TILES = []

assert 'cell_type' in TILE_GDF.columns
assert 'Subtype' in TILE_GDF.columns
assert 'bldgs' in TILE_GDF['Type'].unique()
assert 'prcls' in TILE_GDF['Type'].unique()

for tp in TILE_GDF['cell_type'].unique():
	for st in TILE_GDF[TILE_GDF['cell_type'] == tp]['Subtype'].dropna().unique():
		tile = Tile(
			name=tp,
			buildings=TILE_GDF[(TILE_GDF['Type'] == 'bldgs') & (TILE_GDF['cell_type'] == tp) & (TILE_GDF['Subtype'] == st)],
			parcels=TILE_GDF[(TILE_GDF['Type'] == 'prcls') & (TILE_GDF['cell_type'] == tp) & (TILE_GDF['Subtype'] == st)],
			network=NET_GDF[(NET_GDF['Type'] == tp) & (NET_GDF['Subtype'] == st)],
			trees=TILE_GDF[(TILE_GDF['Type'] == 'trees') & (TILE_GDF['cell_type'] == tp) & (TILE_GDF['Subtype'] == st)],
			bound=TILE_GDF[(TILE_GDF['Type'] == 'bound') & (TILE_GDF['cell_type'] == tp) & (TILE_GDF['Subtype'] == st)],
			subtype=st)
		TILES.append(tile)
		if st in ['Vertical', 'T-Shape', 'L-Middle', 'L-Edge']:
			TILES.append(Tile(tp, tile.buildings, tile.parcels, tile.network, tile.trees, tile.bound, f'{st} Down').flip_vertical())
			TILES.append(Tile(tp, tile.buildings, tile.parcels, tile.network, tile.trees, tile.bound, f'{st} Left').flip_horizontal())
			TILES.append(Tile(tp, tile.buildings, tile.parcels, tile.network, tile.trees, tile.bound, f'{st} Down Left').flip_horizontal().flip_vertical())

"""
TILES = [
	Tile(
		name='Typical Vancouver Single-Family',
		subtype='Standard',
		buildings=TILE_GDF[(TILE_GDF['Type'] == 'bldgs') & (TILE_GDF['cell_type'] == ) & (TILE_GDF['Subtype'] == 'Standard')],
		parcels=TILE_GDF[(TILE_GDF['Type'] == 'prcls') & (TILE_GDF['cell_type']) & (TILE_GDF['Subtype'] == 'Standard')],
		network=NET_GDF[(NET_GDF['Type'] == 'Typical Vancouver Single-Family') & (NET_GDF['Subtype'] == 'Standard')],
		trees=TREE_GDF[(TREE_GDF['cell_type'] == 'Typical_Van_SF') & (TREE_GDF['Subtype'] == 'Standard')]
	),
	Tile(
		name='Typical Vancouver Single-Family (West)',
		subtype='Standard',
		buildings=gpd.read_file(f'{DIR}/Typical_Van_West_SF_bldgs.shp'),
		parcels=gpd.read_file(f'{DIR}/Typical_Van_West_SF_prcls.shp'),
		network=NET_GDF[(NET_GDF['Type'] == 'Typical Vancouver Single-Family (West)') & (NET_GDF['Subtype'] == 'Standard')],
		trees=TREE_GDF[(TREE_GDF['cell_type'] == 'Typical_Van_West_SF') & (TREE_GDF['Subtype'] == 'Standard')]
	),
	Tile(
		name='Treed Large Home',
		subtype='Standard',
		buildings=gpd.read_file(f'{DIR}/Treed_Large_Home_bldgs.shp'),
		parcels=gpd.read_file(f'{DIR}/Treed_Large_Home_prcls.shp'),
		network=NET_GDF[(NET_GDF['Type'] == 'Treed Large Home') & (NET_GDF['Subtype'] == 'Standard')],
		trees=TREE_GDF[(TREE_GDF['cell_type'] == 'Treed_Large_Home') & (TREE_GDF['Subtype'] == 'Standard')]
	),
	Tile(
		name='Coarse Grain',
		subtype='Standard',
		buildings=gpd.read_file(f'{DIR}/Coarse_Grain_bldgs.shp'),
		parcels=gpd.read_file(f'{DIR}/Coarse_Grain_prcls.shp'),
		network=NET_GDF[(NET_GDF['Type'] == 'Coarse Grain') & (NET_GDF['Subtype'] == 'Standard')],
		trees=TREE_GDF[(TREE_GDF['cell_type'] == 'Coarse_Grain') & (TREE_GDF['Subtype'] == 'Standard')]
	),
	Tile(
		name='Open Low Density',
		subtype='Standard',
		buildings=gpd.read_file(f'{DIR}/Open_Low_Density_bldgs.shp'),
		parcels=gpd.read_file(f'{DIR}/Open_Low_Density_prcls.shp'),
		network=NET_GDF[(NET_GDF['Type'] == 'Open Low Density') & (NET_GDF['Subtype'] == 'Standard')],
		trees=TREE_GDF[(TREE_GDF['cell_type'] == 'Open_Low_Density') & (TREE_GDF['Subtype'] == 'Standard')]
	),
    Tile(
        name='Moderate Density',
        subtype='Standard',
        buildings=gpd.read_file(f'{DIR}/Moderate_Density_bldgs.shp'),
        parcels=gpd.read_file(f'{DIR}/Moderate_Density_prcls.shp'),
        network=NET_GDF[(NET_GDF['Type'] == 'Moderate Density') & (NET_GDF['Subtype'] == 'Standard')],
        trees=TREE_GDF[(TREE_GDF['cell_type'] == 'Moderate_Density') & (TREE_GDF['Subtype'] == 'Standard')]
    )]+\
[Tile(
	name='Mid-High Street',
	subtype=subtype,
	buildings=gpd.read_file(f'{DIR}/Mid_High_Street_bldgs.shp'),
	parcels=gpd.read_file(f'{DIR}/Mid_High_Street_prcls.shp'),
	network=NET_GDF[(NET_GDF['Type'] == 'Mid-High Street') & (NET_GDF['Subtype'] == subtype)],
	trees=TREE_GDF[(TREE_GDF['cell_type'] == 'Mid_High_Street') & (TREE_GDF['Subtype'] == subtype)]
) for subtype in VARIATIONS]+\
[Tile(
	name='Dense Nodal',
	subtype=subtype,
	buildings=gpd.read_file(f'{DIR}/Dense_Nodal_bldgs.shp'),
	parcels=gpd.read_file(f'{DIR}/Dense_Nodal_prcls.shp'),
	network=NET_GDF[(NET_GDF['Type'] == 'Dense Nodal') & (NET_GDF['Subtype'] == subtype)],
	trees=TREE_GDF[(TREE_GDF['cell_type'] == 'Dense_Nodal') & (TREE_GDF['Subtype'] == subtype)]
) for subtype in VARIATIONS]
"""