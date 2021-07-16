import os
from morphology.ShapeTools import Shape
import geopandas as gpd
import pandas as pd
from Tiles import TILES
from Tile import Tile
from tqdm import tqdm
from shapely.affinity import scale


class Grid:
	def __init__(self, gdf, tiles, directory='', land_use=None):
		self.gdf = gdf
		self.tiles = tiles
		self.directory = directory
		self.land_use = land_use

	def assert_type_column(self):
		assert 'Type' in self.gdf.columns, KeyError("Type column not found in gdf")

	def assert_high_st_column(self):
		assert 'High St' in self.gdf.columns, KeyError("High St column not found in gdf")

	def assert_land_use_layer(self):
		assert 'elab_lu' in self.land_use.columns, KeyError("elab_lu column not found in land_use layer")

	def assign_subtypes_by_adjency(self, consider_vertices=True):
		subtypes={
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
			# 'L-Middle Left': [{0, 1}, {0, 1, 2, 5}],
			# 'L-Middle': [{1, 2}, {0, 1, 2, 6}, {0, 1, 6}, {0, 1, 2, 4, 5, 6}],
			# 'L-Middle Down Left': [{0, 3}],
			# 'L-Middle Down': [{2, 3}],
			'T-Shape': [
				{1, 2}, {0, 1, 2, 6}, {0, 1, 6}, {0, 1, 2, 4, 5, 6}, {2, 3},
				{0, 1, 2}, {1, 2, 3}, {0, 2, 4}, {0, 1, 2, 6, 7}, {0, 1, 2, 3, 6}, {0, 1, 2, 5, 7},
				{0, 1, 2, 3, 4, 7}, {0, 1, 2, 3, 6, 7}, {0, 1, 2, 3, 4, 5, 7}, {0, 1, 2, 3, 5, 7}, {1, 2, 3, 4, 5, 6},
				{0, 1, 2, 3, 4, 5, 6}, {0, 1, 2, 3, 5, 6, 7}, {0, 1, 2, 3, 4, 5, 6, 7}],
			# 'T-Shape Down': [],
			'T-Shape Left': [
				{0, 1}, {0, 3}, {0, 1, 2, 5}, {0, 1, 3}, {0, 2, 3}, {1, 3, 4}, {0, 1, 2, 4, 5}, {0, 2, 3, 4, 5},
				{0, 2, 4, 6, 7}, {0, 1, 2, 3, 4, 5}, {0, 1, 2, 4, 5, 6, 7}],
		}

		# 'Vertical': [
		# 	{1}, {3},
		# 	{1, 2}, {3, 7}, {1, 3},
		# 	{1, 3, 7}, {1, 3, 6}, {1, 5, 6}, {2, 3, 7}, {1, 5, 6},
		# 	{2, 3, 4, 6}, {2, 3, 6, 7}, {1, 2, 3, 6}, {1, 2, 6, 7}, {2, 3, 5, 7},
		# 	{2, 3, 4, 6, 7}, {1, 2, 3, 4, 7}, {1, 2, 3, 6, 7},
		# 	{1, 2, 3, 5, 6, 7}, {1, 2, 3, 4, 6, 7}],

		"""
		subtypes={
			'Horizontal': [{0}, {2}, {0, 2}, {0, 1, 2, 3}],
			'Vertical Left': [{1, 3, 4}, {0, 3, 4}, {0, 1, 3, 4}, {0, 1, 2, 3, 4}],
			'Vertical Middle': [{1}, {3}, {1, 3}],
			'Vertical': [{1, 3}],
			'L-Middle Left': [{0, 1}],
			'L-Middle': [{1, 2}],
			'L-Middle Down Left': [{0, 3}],
			'L-Middle Down': [{2, 3}],
			'L-Edge Left': [],
			'L-Edge': [],
			'L-Edge Down Left': [],
			'L-Edge Down': [],
			'T-Shape': [{0, 1, 2}, {1, 2, 3}],
			'T-Shape Down': [],
			'T-Shape Left': [{0, 1, 3}, {0, 2, 3}],
		}

		"""

		self.gdf['id'] = self.gdf.index
		gdf = self.gdf.copy()
		type_gdf = self.gdf[self.gdf['High St'] == 1]

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

		gdf['Subtype'] = gdf['Subtype'].fillna('Standard')
		gdf.to_file('del_grid.geojson', driver='GeoJSON')
		return gdf

	def assign_subtypes_by_landuse(self):
		gdf = self.gdf.copy()
		real_land_use_gdf = self.land_use[self.land_use['elab_lu'].isin(['CM', 'IND', 'CV'])].copy()
		real_land_use_gdf['area'] = real_land_use_gdf.area

		for i in tqdm(gdf.index):
			if gdf.loc[i, 'Type'] in ['Open_Low_Density', 'Coarse_Grain']:
				# Get majority of real place land use within the cell
				in_cell_lu = gpd.overlay(real_land_use_gdf, gdf.loc[[i], :].copy())
				in_cell_lu = in_cell_lu[in_cell_lu['elab_lu'].isin(['CV', 'CM', 'IND'])].copy()
				if len(in_cell_lu) > 0:
					larger_use = list(in_cell_lu.groupby('elab_lu').sum().sort_values(by='area', ascending=False).index)[0]
					if larger_use == 'CV': self.gdf.loc[i, 'Subtype'] = 'Standard'
					elif larger_use == 'CM': self.gdf.loc[i, 'Subtype'] = 'Commercial'
					elif larger_use == 'IND': self.gdf.loc[i, 'Subtype'] = 'Industrial'

		return gdf

	def test_assign_subtypes(self):
		assert len(self.gdf.index) == len(self.gdf.index.unique()), IndexError("Duplicated indices found")
		self.assert_type_column()
		self.gdf = self.assign_subtypes_by_adjency()
		if self.land_use is not None:
			self.assert_land_use_layer()
			self.gdf = self.assign_subtypes_by_landuse()
		assert sum(self.gdf['Subtype'].isna()) == 0
		return

	def place_tiles(self):
		"""
		Place tiles in the grid according to their type and export them to directory
		:return:
		"""

		grid = self.gdf.copy()

		bld = gpd.GeoDataFrame()
		pcl = gpd.GeoDataFrame()
		net = gpd.GeoDataFrame()
		veg = gpd.GeoDataFrame()

		# Iterate over tiles
		for tile in tqdm(self.tiles):

			# Filter centroids of subtype
			flt = (grid['Type'] == tile.name) & (grid['Subtype'] == tile.subtype)
			centroids = grid[flt].copy().centroid

			for centroid in centroids:
				tile.all_layers = tile.move_all_layers(centroid)
				bld = pd.concat([bld, tile.all_layers[tile.all_layers['Type'] == 'bldgs']])
				pcl = pd.concat([pcl, tile.all_layers[tile.all_layers['Type'] == 'prcls']])
				net = pd.concat([net, tile.all_layers[tile.all_layers['Type'] == 'ntwkr']])
				veg = pd.concat([veg, tile.all_layers[tile.all_layers['Type'] == 'trees']])

		if len(bld) > 0: bld.to_file(f'{self.directory}/buildings.shp')
		if len(pcl) > 0: pcl.to_file(f'{self.directory}/parcels.shp')
		if len(net) > 0: net.to_file(f'{self.directory}/network.shp')
		if len(veg) > 0: veg.to_file(f'{self.directory}/trees.shp')
		return

	def test_place_tiles(self):
		assert os.path.exists(self.directory), NotADirectoryError(f"{self.directory} not found")
		self.assert_type_column()

		tiles = [tile.name for tile in self.tiles]
		diff = set(self.gdf['Type'].unique()).difference(set(tiles))
		assert len(diff) == 0, KeyError(f"{diff} type not found in Type column ({self.gdf['Type'].unique()})")
		self.place_tiles()
		return self

	def test_grid(self):
		self.test_assign_subtypes()
		self.test_place_tiles()
		return


if __name__ == '__main__':
	TYPES = {
		2: 'Open_Low_Density',
		3: 'Mid_High_Street',
		4: 'Coarse_Grain',
		5: 'Moderate_Density',
		8: 'Treed_Large_Home',
		9: 'Dense_Nodal',
		10: 'Typical_Van_SF',
		11: 'Typical_Van_West_SF'
	}
	grid_gdf = gpd.read_file('/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Sandbox/shp/COV/Data/gmm_grids/fishnet_CoV_gmm_r4_built_broadway.shp')
	grid_gdf['id'] = grid_gdf.index
	grid_gdf['Type'] = grid_gdf['clus_gmm'].replace(TYPES)
	grid_gdf.loc[grid_gdf['Type'].isin(['Mid_High_Street', 'Moderate_Density', 'Dense_Nodal']), 'High St Type'] = 1
	streets = gpd.read_file('/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Sandbox/shp/COV/Data/COV/streets_utm_cov.shp')
	streets['geometry'] = streets.buffer(5)
	grid_gdf.loc[gpd.overlay(grid_gdf, streets[streets['Category'] == 'Arterial'])['id'], 'Arterial'] = 1
	grid_gdf.loc[(grid_gdf['Arterial'] == 1) & (grid_gdf['High St Type'] == 1), 'High St'] = 1
	Grid(grid_gdf, TILES, 'data', land_use=gpd.read_file('/Users/nicholasmartino/Desktop/Landuse2016/Landuse2016.shp')).test_grid()
	# f'/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Sandbox/shp/City-Wide'
