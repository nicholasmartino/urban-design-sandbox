import pandas as pd
from morphology.ShapeTools import Shape
from shapely.affinity import scale


class Tile:
	def __init__(self, name, buildings, parcels, network, trees, bound, subtype=''):
		self.name = name

		for layer, gdf in zip(['buildings', 'parcels', 'network', 'trees'], [buildings, parcels, network, trees]):
			assert len(gdf) > 0, ValueError(f'Empty {layer} GeoDataFrame in {name} tile ({subtype} subtype)')
			assert 'Subtype' in gdf.columns, KeyError(f"Subtype column not found in {layer} layer of {name} tile ({subtype} subtype)")
		assert len(bound) > 0, ValueError(f'Empty boundary GeodataFrame in {name} tile ({subtype} subtype)')

		self.buildings = buildings.copy()
		self.parcels = parcels.copy()
		self.network = network.copy()
		self.trees = trees.copy()
		self.bound = bound.copy()
		self.all_layers = pd.concat([buildings, parcels, network, trees]).reset_index(drop=True)
		self.subtype = subtype

	def move_all_layers(self, reference):
		self.buildings['Type'] = 'bldgs'
		self.parcels['Type'] = 'prcls'
		self.network['Type'] = 'ntwkr'
		self.trees['Type'] = 'trees'

		all_layers = Shape(self.all_layers).move(reference)
		return all_layers

	def flip_horizontal(self):
		self.all_layers = Shape(self.all_layers).flip(-1, 1, self.bound.centroid)
		return self

	def flip_vertical(self):
		self.all_layers = Shape(self.all_layers).flip(1, -1, self.bound.centroid)
		return self
