import numpy as np
import pandas as pd

from Inputs import *
from Network import Streets
from ShapeTools import Shape, Analyst


class Indicators:
	def __init__(self, parcels=None, buildings=None, streets=None, blocks=None):
		parcels['pid'] = parcels.index
		buildings['bid'] = buildings.index
		self.parcels = parcels
		self.buildings = buildings
		self.streets = streets
		self.blocks = blocks
		return

	def get_buildings_floor_area(self):
		gdf = self.buildings.copy()
		gdf['area'] = gdf.area
		gdf['floor_area'] = gdf['area'] * pd.to_numeric(gdf['maxstories'])
		return gdf

	def get_parcel_far(self):
		gdf = self.parcels.copy()
		gdf['floor_area'] = Analyst(gdf, self.buildings.loc[:, ['floor_area', 'geometry']]).spatial_join(operations=['max'])['floor_area_max']
		gdf['area'] = gdf.area
		gdf['parcelarea'] = gdf.area
		gdf['far'] = gdf['floor_area']/gdf['area']
		return gdf

	def get_residential_units(self):

		gdf = self.buildings.copy()
		gdf = gdf.reset_index(drop=True)
		gdf['maxstories'] = gdf['maxstories'].fillna(0)
		gdf['maxstories'] = pd.to_numeric(gdf['maxstories'])

		# Set average unit area by cell type
		sqft_sqm = 10.76391
		gdf.loc[(gdf['cell_type'] == 'Open_Low_Density') & (gdf['laneway'] != '1'), 'unit_area'] = 2257/sqft_sqm
		gdf.loc[(gdf['cell_type'] == 'Mid_High_Street') & (gdf['laneway'] != '1'), 'unit_area'] = 1441/sqft_sqm
		gdf.loc[(gdf['cell_type'] == 'Coarse_Grain') & (gdf['laneway'] != '1'), 'unit_area'] = 2194/sqft_sqm
		gdf.loc[(gdf['cell_type'] == 'Moderate_Density') & (gdf['laneway'] != '1'), 'unit_area'] = 1396/sqft_sqm
		gdf.loc[(gdf['cell_type'] == 'Treed_Large_Home') & (gdf['laneway'] != '1'), 'unit_area'] = 2468/sqft_sqm
		gdf.loc[(gdf['cell_type'] == 'Dense_Nodal') & (gdf['laneway'] != '1'), 'unit_area'] = 1474/sqft_sqm
		gdf.loc[(gdf['cell_type'] == 'Typical_Van_SF') & (gdf['laneway'] != '1'), 'unit_area'] = 1537/sqft_sqm
		gdf.loc[(gdf['cell_type'] == 'Typical_Van_West_SF') & (gdf['laneway'] != '1'), 'unit_area'] = 2058/sqft_sqm
		gdf['unit_area'] = gdf['unit_area'].fillna(0).replace(-np.inf, 0)
		gdf['unit_area'] = gdf['unit_area'] * 1.1
		# Numbers from BC Assessment: Average TOTAL_FINISHED_AREA

		gdf.loc[gdf['LANDUSE'].isin(['CM', 'CV', 'IND', 'OS']), 'res_units'] = 0
		gdf.loc[gdf['LANDUSE'].isin(['MFL', 'MFM', 'MFH', 'SFA']), 'res_units'] = (
				(gdf.area/gdf['unit_area']) * gdf['maxstories']).fillna(0).replace(np.inf, 0).astype(int)
		gdf.loc[(gdf['LANDUSE'] == 'SFD') & (gdf['laneway'] == 0), 'res_units'] = 1
		gdf.loc[gdf['LANDUSE'] == 'MX', 'res_units'] = ((gdf.area/gdf['unit_area']) * (gdf['maxstories'] - 1)).fillna(0).replace(-np.inf, 0).astype(int)
		gdf['res_units'] = gdf['res_units'].fillna(0).astype(int)
		return gdf

	def get_commercial_units(self, unit_area=100):
		gdf = self.buildings.copy()
		gdf['maxstories'] = gdf['maxstories'].fillna(0)

		if sum(gdf['floor_area']) > 0:
			gdf.loc[gdf['LANDUSE'] == 'CM', 'comm_units'] = (gdf['floor_area']/unit_area).astype(int) * gdf['maxstories']
			gdf.loc[gdf['LANDUSE'] == 'MX', 'comm_units'] = (gdf['floor_area']/unit_area).astype(int)
		return gdf

	def get_resident_count(self):
		gdf = self.buildings.copy()
		gdf['res_units'] = pd.to_numeric(gdf['res_units'])
		gdf.loc[gdf['cell_type'] == 'Open_Low_Density', 'res_count'] = 2.098 * gdf.loc[gdf['cell_type'] == 'Open_Low_Density', 'res_units']
		gdf.loc[gdf['cell_type'] == 'Mid_High_Street', 'res_count'] = 1.967 * gdf.loc[gdf['cell_type'] == 'Mid_High_Street', 'res_units']
		gdf.loc[gdf['cell_type'] == 'Coarse_Grain', 'res_count'] = 1.753 * gdf.loc[gdf['cell_type'] == 'Coarse_Grain', 'res_units']
		gdf.loc[gdf['cell_type'] == 'Moderate_Density', 'res_count'] = 1.894 * gdf.loc[gdf['cell_type'] == 'Moderate_Density', 'res_units']
		gdf.loc[gdf['cell_type'] == 'Treed_Large_Home', 'res_count'] = 1.929 * gdf.loc[gdf['cell_type'] == 'Treed_Large_Home', 'res_units']
		gdf.loc[gdf['cell_type'] == 'Dense_Nodal', 'res_count'] = 1.653 * gdf.loc[gdf['cell_type'] == 'Dense_Nodal', 'res_units']
		gdf.loc[gdf['cell_type'] == 'Typical_Van_SF', 'res_count'] = 2.404 * gdf.loc[gdf['cell_type'] == 'Typical_Van_SF', 'res_units']
		gdf.loc[gdf['cell_type'] == 'Typical_Van_West_SF', 'res_count'] = 2.307 * gdf.loc[gdf['cell_type'] == 'Typical_Van_West_SF', 'res_units']
		# Numbers from census canada dissemination areas: Average people per dwelling (Population, 2016 / n_dwellings)
		return gdf

	def remove_buildings_from_open_spaces(self):
		gdf = self.buildings.copy()
		gdf['geometry'] = gdf.centroid
		os_bld = gpd.overlay(gdf, self.parcels[self.parcels['LANDUSE'] == 'OS'])
		return self.buildings[~self.buildings['bid'].isin(os_bld['bid'])]

	def remove_parcels_from_open_spaces(self):
		self.parcels['pid'] = self.parcels.reset_index(drop=True).index
		gdf = self.parcels.copy()
		gdf['geometry'] = gdf.centroid
		try:
			os_pcl = gpd.overlay(gdf, self.parcels[self.parcels['LANDUSE'] == 'OS'].loc[:, ['geometry']])
			return self.parcels[~self.parcels['pid'].isin(os_pcl['pid'])]
		except:
			print("Overlay error, parcels not removed from open spaces")

	def get_area_by_land_use(self):
		gdf = self.parcels.copy()
		gdf['area'] = gdf.area
		df = pd.DataFrame(gdf.groupby('LANDUSE').sum()['area'])
		df['LANDUSE'] = df.index
		df = df.rename(columns={'area': 'Area (m²)', 'LANDUSE': 'Land Use'})
		return df.reset_index(drop=True)

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

	def get_intersection_count(self):
		gdf = self.streets.copy()
		assert 'laneway' in gdf.columns, KeyError("Laneway column not found in streets GeoDataFrame")

		gdf['laneway'] = gdf['laneway'].fillna(0)
		itr = Streets(gdf[gdf['laneway'] == 0].loc[:, ['geometry']]).extract_intersections()
		itr = gpd.overlay(itr, self.parcels.copy(), how="difference")
		return len(itr)

	def get_street_length(self):
		gdf = self.streets.copy()
		assert 'laneway' in gdf.columns, KeyError("Laneway column not found in streets GeoDataFrame")
		gdf['laneway'] = gdf['laneway'].fillna(0)
		gdf = gdf[gdf['laneway'] != "1"]
		parcels = self.parcels.unary_union
		gdf['geometry'] = [geom.difference(parcels) for geom in gdf['geometry']]
		gdf = Streets(gdf).segmentize()
		gdf['centroid'] = [(geom.centroid.x, geom.centroid.y) for geom in gdf['geometry']]
		gdf = gdf.drop_duplicates(['centroid'])
		gdf['length'] = gdf.length
		return gdf

	def get_block_area(self):
		# streets = self.streets.copy()
		# not_lanes = self.streets[(self.streets['laneway'] == 0) | (self.streets['laneway'].isna())].copy()
		# not_lanes['geometry'] = not_lanes.buffer(6)
		#
		# not_lanes_os = gpd.overlay(not_lanes, self.parcels.loc[self.parcels['LANDUSE'] == 'OS', ['geometry']].copy(), how='difference')
		# blocks_gdf = gpd.overlay(gpd.GeoDataFrame({'geometry': [self.parcels.unary_union.convex_hull]}, crs=self.parcels.crs), not_lanes_os, how='difference')
		# blocks_gdf = Shape(blocks_gdf).divorce()

		# # Buffer laneways
		# lanes = streets[streets['laneway'] == '1'].copy()
		# lanes['geometry'] = lanes.buffer(6)
		#
		# # Remove lanes that intersect with open spaces
		# parcels = self.parcels.copy()
		# parcels['geometry'] = parcels.buffer(5)
		# lanes = gpd.overlay(lanes.loc[:, ['geometry']], parcels[parcels['LANDUSE'] == 'OS'], how="difference")
		# lanes = gpd.overlay(lanes.loc[:, ['geometry']], self.blocks)

		# Merge parcels and buffered lanes to form blocks
		parcels = self.parcels.copy()
		parcels.loc[parcels['LANDUSE'] == 'OS', 'geometry'] = [geom.buffer(5) for geom in parcels.loc[parcels['LANDUSE'] == 'OS', 'geometry']]
		blocks_gdf = gpd.overlay(self.blocks.copy(), parcels[parcels['LANDUSE'] == 'OS'], how="difference")
		blocks_gdf = Shape(pd.concat([self.parcels[self.parcels['LANDUSE'] == 'OS'], blocks_gdf])).dissolve()
		blocks_gdf = Shape(blocks_gdf).divorce()
		blocks_gdf['area'] = blocks_gdf.area
		return blocks_gdf

	def get_height_from_storeis(self):
		gdf = self.buildings.copy()
		gdf['height'] = gdf['maxstories'] * 3
		return gdf

	def join_land_use_from_parcels(self):
		gdf = self.buildings.copy().reset_index(drop=True)
		gdf['LANDUSE'] = Analyst(gdf, self.parcels.copy().loc[:, ['LANDUSE', 'geometry']]).spatial_join()['LANDUSE']
		return gdf

	def test_indicators(self):
		bld_cols = ['LANDUSE', 'maxstories', 'laneway']
		for col in bld_cols:
			assert col in self.buildings.columns, KeyError(f'{col} column not found in buildings layer')
		self.buildings['maxstories'] = self.buildings['maxstories'].fillna(0)
		self.buildings['maxstories'] = pd.to_numeric(self.buildings['maxstories'])
		assert 'LANDUSE' in self.parcels.columns, KeyError('LANDUSE column not found in parcels layer')
		assert 'pid' in self.parcels.columns, KeyError('pid column not found in parcels layer')

		print("Extracting indicators")
		self.buildings = self.get_buildings_floor_area()
		self.parcels = self.get_parcel_far()
		self.buildings = self.get_residential_units()
		self.buildings = self.get_commercial_units()
		self.buildings = self.get_resident_count()
		self.buildings = self.get_height_from_storeis()
		print("Indicators successfully extracted")

		if len(self.parcels[self.parcels['LANDUSE'] == 'OS']) > 0:
			self.buildings = self.remove_buildings_from_open_spaces()
		# self.parcels = self.remove_parcels_from_open_spaces()
		self.buildings = self.join_land_use_from_parcels()
		self.get_area_by_land_use()

		return


class Scenario:
	def __init__(self, parcels, buildings, trees, real_parks, real_trees, name=''):
		self.parcels = parcels
		self.trees = trees
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
		prcls.loc[prcls['LANDUSE'] == 'OS', ['geometry']] = prcls.loc[prcls['LANDUSE'] == 'OS', ['geometry']].buffer(12, join_style=2)
		if len(prcls[prcls['LANDUSE'] == 'OS']) > 0:
			dissolved = Shape(prcls[prcls['LANDUSE'] == 'OS']).dissolve()
			dissolved['LANDUSE'] = 'OS'
			dissolved['Type'] = 'prcls'
			prcls = pd.concat([prcls[prcls['LANDUSE'] != 'OS'], dissolved])
			prcls.loc[prcls['LANDUSE'] == 'OS', ['geometry']] = prcls.loc[prcls['LANDUSE'] == 'OS', ['geometry']].buffer(-16, join_style=2)
			prcls = prcls.dropna(subset=['geometry']).reset_index(drop=True)
		return prcls

	def extract_trees(self, directory=''):
		trees = self.real_trees

		assert 'crown_dm' in trees.columns, KeyError("'Crown_dm' column not found in trees GeoDataFrame")
		assert 'LANDUSE' in self.parcels.columns, KeyError("'LANDUSE' column not found in parcels GeoDataFrame")

		# Get real place trees
		trees['real_tree_id'] = trees.index
		trees['Crown_dm'] = trees['crown_dm']
		real_trees_copy = trees.copy()
		real_trees_copy['geometry'] = trees.centroid.buffer(1)
		overlay = gpd.overlay(real_trees_copy.loc[:, ['real_tree_id', 'geometry']], self.parcels[self.parcels['LANDUSE'] == 'OS'])
		sb_trees = self.trees.copy()
		sb_trees = pd.concat([sb_trees, trees[trees['real_tree_id'].isin(list(overlay['real_tree_id']))]])
		sb_trees['Type'] = 'trees'
		return sb_trees
