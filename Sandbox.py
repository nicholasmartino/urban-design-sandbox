import pandas as pd
import os
import geopandas as gpd
from morphology.ShapeTools import Shape, Analyst


class Indicators:
	def __init__(self, parcels, buildings):
		parcels['pid'] = parcels.index
		buildings['bid'] = buildings.index
		self.parcels = parcels
		self.buildings = buildings
		return

	def get_buildings_floor_area(self):
		gdf = self.buildings.copy()
		gdf['area'] = gdf.area
		gdf['floor_area'] = gdf['area'] * gdf['maxstories']
		return gdf

	def get_parcel_far(self):
		gdf = self.parcels.copy()
		gdf['floor_area'] = Analyst(gdf, self.buildings.loc[:, ['floor_area', 'geometry']]).spatial_join(operations=['max'])['floor_area_max']
		gdf['area'] = gdf.area
		gdf['parcelarea'] = gdf.area
		gdf['far'] = gdf['floor_area']/gdf['area']
		return gdf

	def get_residential_units(self, unit_area=270):
		gdf = self.buildings.copy()
		gdf.loc[gdf['LANDUSE'].isin(['CM', 'CV', 'IND', 'OS']), 'res_units'] = 0
		gdf.loc[gdf['LANDUSE'].isin(['MFL', 'MFM', 'MFH', 'SFA']), 'res_units'] = ((gdf['area']/unit_area) * gdf['maxstories']).astype(int)
		gdf.loc[(gdf['LANDUSE'] == 'SFD') & (gdf['laneway'] == 0), 'res_units'] = 1
		gdf.loc[gdf['LANDUSE'] == 'MX', 'res_units'] = ((gdf['area']/unit_area) * (gdf['maxstories'] - 1)).astype(int)
		gdf['res_units'] = gdf['res_units'].fillna(0)
		return gdf

	def get_commercial_units(self, unit_area=100):
		gdf = self.buildings.copy()
		gdf.loc[gdf['LANDUSE'] == 'CM', 'comm_units'] = (gdf['floor_area']/unit_area).astype(int) * gdf['maxstories']
		gdf.loc[gdf['LANDUSE'] == 'MX', 'comm_units'] = (gdf['floor_area']/unit_area).astype(int)
		return gdf

	def get_resident_count(self):
		gdf = self.buildings.copy()
		gdf['res_count'] = gdf['res_units'].astype(int) * 3
		return gdf

	def remove_buildings_from_open_spaces(self):
		gdf = self.buildings.copy()
		gdf['geometry'] = gdf.centroid
		os_bld = gpd.overlay(gdf, self.parcels[self.parcels['LANDUSE'] == 'OS'])
		return self.buildings[~self.buildings['bid'].isin(os_bld['bid'])]

	def get_area_by_land_use(self):
		gdf = self.parcels.copy()
		return pd.DataFrame(gdf.groupby('LANDUSE').sum()['area'])

	def get_floor_area_by_land_use(self):
		gdf = self.buildings.copy()
		return pd.DataFrame(gdf.groupby('LANDUSE').sum()['floor_area'])

	def get_n_units_by_land_use(self):
		gdf = self.buildings.copy()
		gdf['res_units'] = gdf['res_units'].astype(int)
		return pd.DataFrame(gdf.groupby('LANDUSE').sum()['res_units'])

	def get_total_population(self):
		gdf = self.buildings.copy()
		return int(gdf['res_count'].sum())

	def get_total_area(self):
		gdf = self.parcels.copy()
		return round(gdf.unary_union.convex_hull.area, 2)

	def join_land_use_from_parcels(self):
		gdf = self.buildings.copy()
		gdf['LANDUSE'] = gpd.sjoin(gdf, self.parcels.copy().loc[:, ['LANDUSE', 'geometry']])['LANDUSE_right']
		return gdf

	def test_indicators(self):
		bld_cols = ['LANDUSE', 'maxstories', 'laneway']
		for col in bld_cols:
			assert col in self.buildings.columns, KeyError(f'{col} column not found in buildings layer')
		self.buildings['maxstories'] = pd.to_numeric(self.buildings['maxstories'])
		assert 'LANDUSE' in self.parcels.columns, KeyError('LANDUSE column not found in parcels layer')
		assert 'pid' in self.parcels.columns, KeyError('pid column not found in parcels layer')

		self.buildings = self.get_buildings_floor_area()
		self.parcels = self.get_parcel_far()
		self.buildings = self.get_residential_units()
		self.buildings = self.get_commercial_units()
		self.buildings = self.get_resident_count()

		assert len(self.parcels[self.parcels['LANDUSE'] == 'OS']) > 0, AssertionError('Parcel layer has no open spaces (OS) under LANDUSE column')
		self.buildings = self.remove_buildings_from_open_spaces()
		self.buildings = self.join_land_use_from_parcels()
		self.get_area_by_land_use()

		return


class Scenario:
	def __init__(self, parcels, buildings, real_parks, real_trees, name=''):
		self.parcels = parcels
		self.buildings = buildings
		self.real_parks = real_parks
		self.real_trees = real_trees
		self.name = name
		return

	def extract_parks(self):
		prcls = self.parcels.copy()

		assert 'LANDUSE' in prcls.columns, KeyError("'LANDUSE' column not found in parcels GeoDataFrame")

		# Extract open spaces from real place
		prcls['pid'] = prcls.index
		pcl = prcls.copy()
		pcl['geometry'] = pcl.centroid
		inters_parcels = gpd.overlay(pcl.loc[:, ['pid', 'geometry']], self.real_parks)
		prcls.loc[prcls['pid'].isin(inters_parcels['pid']), 'LANDUSE'] = 'OS'
		prcls.loc[prcls['LANDUSE'] == 'OS', ['geometry']] = prcls.loc[prcls['LANDUSE'] == 'OS', ['geometry']].buffer(5, join_style=2)
		dissolved = Shape(prcls[prcls['LANDUSE'] == 'OS']).dissolve()
		dissolved['LANDUSE'] = 'OS'
		dissolved['Type'] = 'prcls'
		prcls = pd.concat([prcls[prcls['LANDUSE'] != 'OS'], dissolved])
		prcls.loc[prcls['LANDUSE'] == 'OS', ['geometry']] = prcls.loc[prcls['LANDUSE'] == 'OS', ['geometry']].buffer(-6, join_style=2)
		prcls = prcls.dropna(subset=['geometry']).reset_index(drop=True)
		return prcls

	def extract_trees(self):
		trees = self.real_trees

		assert 'crown_dm' in trees.columns, KeyError("'Crown_dm' column not found in trees GeoDataFrame")
		assert 'LANDUSE' in self.parcels.columns, KeyError("'LANDUSE' column not found in parcels GeoDataFrame")

		assert os.path.exists(f'data/sandboxes/{self.name}trees.shp'), FileNotFoundError(f"'data/sandboxes/{self.name}trees.shp' file not found")

		# Get real place trees
		trees['real_tree_id'] = trees.index
		trees['Crown_dm'] = trees['crown_dm']
		real_trees_copy = trees.copy()
		real_trees_copy['geometry'] = trees.centroid.buffer(1)
		overlay = gpd.overlay(real_trees_copy.loc[:, ['real_tree_id', 'geometry']], self.parcels[self.parcels['LANDUSE'] == 'OS'])
		sb_trees = gpd.read_file(f'data/sandboxes/{self.name}trees.shp')
		sb_trees = pd.concat([sb_trees, trees[trees['real_tree_id'].isin(list(overlay['real_tree_id']))]])
		sb_trees['Type'] = 'trees'
		return sb_trees


if __name__ == '__main__':
	PREFIXES = ['broadway_baseline.geojson', 'broadway_e1.geojson']

	# Get open spaces from CoV open data
	PARKS = gpd.read_file('https://opendata.vancouver.ca/explore/dataset/parks-polygon-representation/download/'
	                      '?format=geojson&timezone=America/Los_Angeles&lang=en&').to_crs(26910)
	REAL_TREES = gpd.read_file('/Volumes/Samsung_T5/Databases/Vancouver, British Columbia.gpkg',
	                           layer='nature_tree_canopy')

	for prefix in PREFIXES:
		prefix = f"{prefix.split('.')[0]}_"
		prcls = gpd.read_file(f"data/sandboxes/{prefix}parcels.shp")
		bldgs = gpd.read_file(f"data/sandboxes/{prefix}buildings.shp")

		scn = Scenario(parcels=prcls, buildings=bldgs, real_parks=PARKS, real_trees=REAL_TREES, name=prefix)
		scn.parcels = scn.extract_parks()
		sb_trees = scn.extract_trees()

		ind = Indicators(parcels=scn.parcels, buildings=bldgs)
		ind.test_indicators()
		ind.parcels.to_file(f"data/sandboxes/{prefix}parcels_indicator.shp")
		ind.buildings.to_file(f"data/sandboxes/{prefix}buildings_indicator.shp")
		ind.get_area_by_land_use().to_csv(f'data/sandboxes/{prefix}land_use_area.csv')
		ind.get_floor_area_by_land_use().to_csv(f'data/sandboxes/{prefix}land_use_floor_area.csv')
		ind.get_n_units_by_land_use().to_csv(f'data/sandboxes/{prefix}land_use_n_units.csv')

		all_tiles = pd.concat([ind.parcels, ind.buildings, sb_trees])
		all_tiles.to_file(f'data/sandboxes/{prefix}all_tiles.geojson', driver='GeoJSON')
		all_tiles.to_crs(4326).to_file(f'data/sandboxes/{prefix}all_tiles_4326.geojson', driver='GeoJSON')

		print(f"Total population: {ind.get_total_population()}")
		print(f"Total area: {ind.get_total_area()}")
