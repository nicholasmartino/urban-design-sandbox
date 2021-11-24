import pandas as pd
import sys
from Grid import Grid
from Inputs import *
from Sandbox import Scenario, Indicators
from Tiles import TILES


TYPES = {
	2: 'Open_Low_Density',
	3: 'Mid_High_Street',
	4: 'Coarse_Grain',
	5: 'Moderate_Density',
	8: 'Treed_Large_Home',
	9: 'Dense_Nodal',
	10: 'Typical_Van_SF',
	11: 'Typical_Van_West_SF',
	12: 'Green_Open_Space'
}

for file in GRID_FILES:
	grid_file = f'{GRID_DIR}/{file}'
	grid_gdf = gpd.read_file(grid_file)
	grid_gdf['id'] = grid_gdf.index
	grid_gdf['Type'] = grid_gdf['clus_gmm'].replace(TYPES)
	grid_gdf.loc[grid_gdf['Type'].isin(['Mid_High_Street', 'Moderate_Density', 'Dense_Nodal']), 'High St Type'] = 1
	STREETS['geometry'] = STREETS.buffer(5)
	grid_gdf.loc[gpd.overlay(grid_gdf, STREETS[STREETS['Category'] == 'Arterial'])['id'], 'Arterial'] = 1
	grid_gdf.loc[(grid_gdf['Arterial'] == 1) & (grid_gdf['High St Type'] == 1), 'High St'] = 1
	Grid(gdf=grid_gdf, tiles=TILES, directory=OUT_DIR, prefix=f"{file.split('.')[0]}_",
	     land_use=gpd.read_file('data/mvan/Landuse2016/Landuse2016.shp'),
	     diagonal_gdf=gpd.read_file('data/diagonal_tiles.geojson')).test_grid()

# Get open spaces from CoV open data
if sys.platform == 'win32':
	GPKG = 'G:/My Drive/Databases/Vancouver, British Columbia.gpkg'
else:
	GPKG = '/Volumes/Samsung_T5/Databases/Vancouver, British Columbia.gpkg'
PARKS = gpd.read_file('https://opendata.vancouver.ca/explore/dataset/parks-polygon-representation/download/'
                      '?format=geojson&timezone=America/Los_Angeles&lang=en&').to_crs(26910)
REAL_TREES = gpd.read_file(GPKG, layer='nature_tree_canopy')

for prefix in GRID_FILES:
	prefix = f"{prefix.split('.')[0]}_"
	print(f"\n {prefix}")
	prcls = gpd.read_file(f"{OUT_DIR}/{prefix}parcels.shp")
	bldgs = gpd.read_file(f"{OUT_DIR}/{prefix}buildings.shp")
	strts = gpd.read_file(f"{OUT_DIR}/{prefix}network.shp")
	blcks = gpd.read_file(f"{OUT_DIR}/{prefix}blocks.shp")

	scn = Scenario(parcels=prcls, buildings=bldgs, real_parks=PARKS, real_trees=REAL_TREES, name=prefix)
	scn.parcels = scn.extract_parks()
	sb_trees = scn.extract_trees(directory=OUT_DIR)

	ind = Indicators(parcels=scn.parcels, buildings=bldgs, streets=strts, blocks=blcks)
	ind.test_indicators()
	ind.parcels.to_file(f"{OUT_DIR}/{prefix}parcels_indicator.shp")
	ind.buildings.to_file(f"{OUT_DIR}/{prefix}buildings_indicator.shp")
	ind.get_area_by_land_use().to_csv(f'{OUT_DIR}/{prefix}land_use_area.csv')
	ind.get_floor_area_by_land_use().to_csv(f'{OUT_DIR}/{prefix}land_use_floor_area.csv')
	ind.get_n_units_by_land_use().to_csv(f'{OUT_DIR}/{prefix}land_use_n_units.csv')
	print(f"Intersection count: {ind.get_intersection_count()}")
	# stt = ind.get_street_length()
	# stt.loc[:, ["length", "geometry"]].to_file(f'{OUT_DIR}/{prefix}network_indicator.shp')
	# stt.loc[:, ["length"]].to_csv(f'{OUT_DIR}/{prefix}street_length.csv')
	# blocks = ind.get_block_area()
	# blocks.to_csv(f'{OUT_DIR}/{prefix}block_area.csv')
	# blocks.to_file(f'{OUT_DIR}/{prefix}blocks_indicator.shp')

	all_tiles = pd.concat([ind.parcels, ind.buildings, sb_trees])
	all_tiles.to_file(f'{OUT_DIR}/{prefix}all_tiles.geojson', driver='GeoJSON')
	# all_tiles.to_crs(4326).to_file(f'{OUT_DIR}/{prefix}all_tiles_4326.geojson', driver='GeoJSON')

	print(f"Total population: {ind.get_total_population()}")
	print(f"Total area: {ind.get_total_area()}")
