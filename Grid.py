import os

import geopandas as gpd
import pandas as pd
from Tiles import TILES, DIR
from tqdm import tqdm


class Grid:
	def __init__(self, gdf, tiles, directory=''):
		self.gdf = gdf
		self.tiles = tiles
		self.directory = directory

	def place_tiles(self):
		"""
		Place tiles in the grid according to their type and export them to directory
		:param directory: Export directory
		:return:
		"""

		grid = self.gdf.copy()

		# Iterate over tiles
		for tile in tqdm(self.tiles):
			bld = gpd.GeoDataFrame()
			pcl = gpd.GeoDataFrame()
			net = gpd.GeoDataFrame()

			# Filter grid for type
			centroids = grid[grid['Type'] == tile.name].centroid

			# Iterate over centroids
			for point in centroids:
				tile.move_all_layers(point)
				bld = pd.concat([bld, tile.buildings])
				pcl = pd.concat([pcl, tile.parcels])
				net = pd.concat([net, tile.network])

			bld.to_file(f'{self.directory}/{tile.name}_bld.shp')
			pcl.to_file(f'{self.directory}/{tile.name}_pcl.shp')
			if len(net) > 0: net.to_file(f'{self.directory}/{tile.name}_net.shp')
		return

	def test_place_tiles(self):
		assert os.path.exists(self.directory), NotADirectoryError(f"{self.directory} not found")
		assert 'Type' in self.gdf.columns, KeyError("Type column not found in gdf")
		for tile in self.tiles:
			assert tile.name in self.gdf['Type'].unique(), KeyError(f"{tile.name} type not found in Type column")
		self.place_tiles()


if __name__ == '__main__':
	TYPES = {
		2: 'Open Low Density',
		3: 'Mid-High Street',
		4: 'Coarse Grain',
		5: 'Moderate Density',
		8: 'Treed Large Home',
		9: 'Dense Nodal',
		10: 'Typical Vancouver Single-Family',
		11: 'Typical Vancouver Single-Family (West)'
	}
	grid_gdf = gpd.read_file('/Volumes/ELabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Sandbox/shp/COV/Data/gmm_grids/fishnet_CoV_gmm_r4_built.shp')
	grid_gdf['Type'] = grid_gdf['clus_gmm'].replace(TYPES)
	Grid(grid_gdf, TILES, f'/Volumes/ELabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Sandbox/shp/City-Wide').test_place_tiles()
