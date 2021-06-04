import pandas as pd
from Morphology.ShapeTools import Shape


class Tile:
	def __init__(self, name, buildings, parcels, network, trees, subtype=''):
		self.name = name
		self.subtype = subtype

		for layer, gdf in zip(['buildings', 'parcels', 'network', 'trees'], [buildings, parcels, network, trees]):
			assert 'Subtype' in gdf.columns, KeyError(f"Subtype column not found in {layer} layer of {name} tile")
			assert subtype in list(gdf['Subtype']), KeyError(f"{subtype} subtype not found in {layer} layer of {name} tile")

		self.buildings = buildings[buildings['Subtype'] == subtype]
		self.parcels = parcels[parcels['Subtype'] == subtype]
		self.network = network[network['Subtype'] == subtype]
		self.trees = trees[trees['Subtype'] == subtype]

	def move_all_layers(self, reference):
		self.buildings['Type'] = 'bldgs'
		self.parcels['Type'] = 'prcls'
		self.network['Type'] = 'ntwkr'
		self.trees['Type'] = 'trees'

		all_layers = Shape(pd.concat([self.buildings, self.parcels, self.network, self.trees]), reset_index=True).move(reference)
		return all_layers
