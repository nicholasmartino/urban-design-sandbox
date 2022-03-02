import gc

import pandas as pd
from tqdm import tqdm

from models.ShapeTools import Shape, Analyst
from store import *


class Grid:
	def __init__(self, gdf, tiles, directory=os.getcwd(), prefix='', land_use=None, diagonal_gdf=None):
		self.gdf = gdf
		self.tiles = tiles
		self.directory = directory
		self.land_use = land_use
		self.diagonal_gdf = diagonal_gdf
		self.prefix = prefix

		assert os.path.exists(directory), NotADirectoryError(f'{directory}')

	def assert_type_column(self):
		assert 'Type' in self.gdf.columns, KeyError("Type column not found in gdf")

	def assert_subtype_column(self):
		assert 'Subtype' in self.gdf.columns, KeyError("Subtype column not found in gdf")

	def assert_high_st_column(self):
		assert 'High St' in self.gdf.columns, KeyError("High St column not found in gdf")

	def assert_land_use_layer(self):
		assert 'LANDUSE' in self.land_use.columns, KeyError("LANDUSE column not found in land_use layer")

	def assign_subtypes_by_adjency(self, consider_vertices=True):
		subtypes = {
			'Horizontal': [
				{0}, {2}, {5}, {6}, {7},
				{0, 2}, {4, 7},
				{0, 2, 7}, {0, 2, 5}, {0, 2, 6},
				{0, 2, 4, 7}, {0, 2, 4, 5}, {0, 1, 2, 3}, {0, 2, 4, 7}, {0, 2, 5, 6}, {0, 2, 6, 7},
				{0, 2, 4, 5, 6}],
			'Diagonal': [],
			'Diagonal Cross': [],
			'Standard Left': [
				{4},
				{1, 5},
				{1, 3, 5}, {0, 3, 4}, {0, 1, 5}, {3, 4}, {3, 4, 6},
				{0, 3, 4, 7}, {0, 1, 3, 5}, {0, 1, 4, 5}, {0, 1, 3, 4}, {0, 1, 4, 7}, {0, 3, 4, 5},
				{0, 3, 4, 5, 7}, {0, 1, 3, 4, 5}, {0, 1, 2, 3, 4}, {1, 3, 4, 5, 6},
				{0, 1, 3, 4, 5, 6}, {0, 1, 3, 4, 5, 7}],
			'T-Shape': [
				{1, 2}, {0, 1, 2, 6}, {0, 1, 6}, {0, 1, 2, 4, 5, 6}, {2, 3},
				{0, 1, 2}, {1, 2, 3}, {0, 2, 4}, {0, 1, 2, 6, 7}, {0, 1, 2, 3, 6}, {0, 1, 2, 5, 7},
				{0, 1, 2, 3, 4, 7}, {0, 1, 2, 3, 6, 7}, {0, 1, 2, 3, 4, 5, 7}, {0, 1, 2, 3, 5, 7}, {1, 2, 3, 4, 5, 6},
				{0, 1, 2, 3, 4, 5, 6}, {0, 1, 2, 3, 5, 6, 7}, {0, 1, 2, 3, 4, 5, 6, 7}],
			'T-Shape Left': [
				{0, 1}, {0, 3}, {0, 1, 2, 5}, {0, 1, 3}, {0, 2, 3}, {1, 3, 4}, {0, 1, 2, 4, 5}, {0, 2, 3, 4, 5},
				{0, 2, 4, 6, 7}, {0, 1, 2, 3, 4, 5}, {0, 1, 2, 4, 5, 6, 7}],
		}

		self.gdf['id'] = self.gdf.index
		gdf = self.gdf.copy()
		type_gdf = self.gdf[self.gdf['High St'] == 1]
		if len(type_gdf) > 0:
			for i in tqdm(type_gdf.index):

				# Filter grid and buffer cells to overlay with adjacents
				gdf.at[i, 'geometry'] = gdf.loc[i, 'geometry'].buffer(3, join_style=2)

				# Explode cell into segments and centroids
				edges = Shape(Shape(gdf.loc[[i], ['geometry']]).explode()).divorce()
				edges['geometry'] = edges.centroid
				if consider_vertices:
					vertices = Shape(gdf.loc[[i], ['geometry']]).extract_vertices().drop_duplicates(subset=['geometry'])
					vertices['sid'] = range(4, 8)
					points = pd.concat([edges, vertices]).reset_index(drop=True)
				else:
					points = edges.copy()
				points['geometry'] = points.buffer(0.1)

				# Intersect edges and vertices with adjacent geometries
				inters = gpd.overlay(type_gdf, points.loc[:, ['geometry', 'sid']])
				adj_cells = type_gdf[type_gdf['id'].isin(inters['id'])].copy()
				adj_cells['sid'] = list(inters['sid'])

				# Assign subtype based on adjacency patterns
				if len(adj_cells) == 0:
					gdf.loc[i, 'Subtype'] = 'Standard'
				else:
					# Classify subtype according to adjacent cells
					for sub, sets in subtypes.items():
						for val in sets:
							if set(adj_cells['sid']) == val:
								gdf.loc[i, 'Subtype'] = sub

				gdf.at[i, 'geometry'] = gdf.loc[i, 'geometry'].buffer(-3, join_style=2)

			gdf['Subtype'] = gdf['Subtype'].fillna('Standard')
			return gdf
		else:
			print("No High St type found")
			return gdf

	def get_real_place_land_use(self, land_uses):
		real_land_use_gdf = self.land_use[self.land_use['LANDUSE'].isin(land_uses)].copy()
		real_land_use_gdf['area'] = real_land_use_gdf.area
		return real_land_use_gdf

	def get_majority_land_use_in_cell(self, real_land_use_gdf, cell_i, land_uses):
		gdf = self.gdf.copy()
		# Get majority of real place land use within the cell
		in_cell_lu = gpd.overlay(real_land_use_gdf, gdf.loc[[cell_i], :].copy())
		in_cell_lu = in_cell_lu[in_cell_lu['LANDUSE'].isin(land_uses)].copy()
		in_cell_lu['area'] = in_cell_lu.area
		if len(in_cell_lu) > 0:
			return list(in_cell_lu.groupby('LANDUSE').sum().sort_values(by='area', ascending=False).index)[0]

	def assign_subtypes_by_landuse(self):
		gdf = self.gdf.copy()
		real_land_use_gdf_1 = self.get_real_place_land_use(['CV', 'CM', 'IND'])
		real_land_use_gdf_2 = self.get_real_place_land_use(['SFD', 'MFL'])

		for i in tqdm(gdf.index):
			if gdf.loc[i, 'Type'] in ['Open_Low_Density']:
				larger_use = self.get_majority_land_use_in_cell(real_land_use_gdf_1, i, ['CV', 'CM', 'IND'])
				if larger_use == 'CV': gdf.loc[i, 'Subtype'] = 'Standard'
				elif larger_use == 'CM': gdf.loc[i, 'Subtype'] = 'Commercial'
				elif larger_use == 'IND': gdf.loc[i, 'Subtype'] = 'Industrial'

			elif gdf.loc[i, 'Type'] in ['Treed_Large_Home']:
				larger_use = self.get_majority_land_use_in_cell(real_land_use_gdf_2, i, ['SFD', 'MFL'])
				if larger_use == 'SFD': gdf.loc[i, 'Subtype'] = 'Standard'
				elif larger_use == 'MFL': gdf.loc[i, 'Subtype'] = 'Multi-Family'

		if self.diagonal_gdf is not None:
			self.diagonal_gdf['geometry'] = self.diagonal_gdf.buffer(1)
			ids = gpd.overlay(gdf.loc[:, ['id', 'geometry']], self.diagonal_gdf.loc[:, ['geometry']])['id']
			gdf.loc[list(ids), 'Type'] = 'Mid_High_Street'
			gdf.loc[list(ids), 'Subtype'] = 'Diagonal'
		return gdf

	def place_tiles(self, export=False):
		"""
		Place tiles in the grid according to their type and export them to directory
		:return:
		"""

		grid = self.gdf.copy()

		bld = gpd.GeoDataFrame()
		pcl = gpd.GeoDataFrame()
		net = gpd.GeoDataFrame()
		veg = gpd.GeoDataFrame()
		blk = gpd.GeoDataFrame()

		# Iterate over tiles
		all_layers = gpd.GeoDataFrame()
		for tile in tqdm(self.tiles):

			# Filter centroids of subtype
			flt = (grid['Type'] == tile.name) & (grid['Subtype'] == tile.subtype)
			centroids = grid[flt].copy().centroid

			for centroid in centroids:
				tile.all_layers = tile.move_all_layers(centroid)
				tile.all_layers['id_grid'] = list(Analyst(tile.all_layers, grid.loc[:, ['id', 'geometry']]).spatial_join(operations=['max'])['id_max'])
				bld = pd.concat([bld, tile.all_layers[tile.all_layers['Type'] == 'bldgs']])
				pcl = pd.concat([pcl, tile.all_layers[tile.all_layers['Type'] == 'prcls']])
				net = pd.concat([net, tile.all_layers[tile.all_layers['Type'] == 'ntwrk']])
				veg = pd.concat([veg, tile.all_layers[tile.all_layers['Type'] == 'trees']])
				blk = pd.concat([blk, tile.all_layers[tile.all_layers['Type'] == 'block']])
				all_layers = pd.concat([all_layers, tile.all_layers])
				gc.collect()

		if export:
			if len(bld) > 0:
				bld.to_file(f'{self.directory}/{self.prefix}buildings.shp')
			if len(pcl) > 0:
				pcl.to_file(f'{self.directory}/{self.prefix}parcels.shp')
			if len(net) > 0:
				grid_bnd = grid.copy()
				grid_bnd['geometry'] = grid_bnd['geometry'].boundary
				grid_bnd = Shape(grid_bnd).explode()
				net = pd.concat([net, grid_bnd.loc[:, ['geometry']]])
				net = net.drop_duplicates('geometry')
				net.to_file(f'{self.directory}/{self.prefix}network.shp')
			if len(veg) > 0:
				veg.to_file(f'{self.directory}/{self.prefix}trees.shp')
			if len(blk) > 0:
				blk.to_file(f'{self.directory}/{self.prefix}blocks.shp')
		# all_layers = all_layers.reset_index(drop=True).drop('id', axis=1).dropna(how='all')
		# all_layers.to_file(f'{self.directory}/{self.prefix}all_tiles.geojson', driver='GeoJSON')
		# all_layers.to_crs(4326).to_file(f'{self.directory}/{self.prefix}all_tiles4326.geojson', driver='GeoJSON')
		return all_layers

	def test_assign_subtypes(self):
		assert len(self.gdf.index) == len(self.gdf.index.unique()), IndexError("Duplicated indices found")
		self.assert_type_column()
		self.assert_high_st_column()
		self.gdf = self.assign_subtypes_by_adjency()
		if self.land_use is not None:
			self.assert_land_use_layer()
			self.gdf = self.assign_subtypes_by_landuse()
		assert sum(self.gdf['Subtype'].isna()) == 0
		return self.gdf

	def test_place_tiles(self):
		assert os.path.exists(self.directory), NotADirectoryError(f"{self.directory} not found")
		self.assert_type_column()

		tiles = [tile.name for tile in self.tiles]
		diff = set(self.gdf['Type'].unique()).difference(set(tiles))
		assert len(diff) == 0, KeyError(f"{diff} type not found in Type column ({self.gdf['Type'].unique()})")
		tiles = self.place_tiles()
		return tiles

	def test_grid(self):
		self.test_assign_subtypes()
		all_layers = self.test_place_tiles().gdf
		return all_layers
