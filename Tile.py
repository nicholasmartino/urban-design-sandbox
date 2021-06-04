from Morphology.ShapeTools import Shape


class Tile:
	def __init__(self, name, buildings, parcels, network, subtype=''):
		self.name = name
		self.subtype = subtype

		assert 'Subtype' in buildings.columns, KeyError(f"Subtype column not found in buildings layer of {name} tile")
		assert 'Subtype' in parcels.columns, KeyError(f"Subtype column not found in parcels layer of {name} tile")
		assert 'Subtype' in network.columns, KeyError(f"Subtype column not found in network layer of {name} tile")
		assert subtype in list(buildings['Subtype']), KeyError(f"{subtype} subtype not found in buildings layer of {name} tile")
		assert subtype in list(parcels['Subtype']), KeyError(f"{subtype} subtype not found in parcels layer of {name} tile")
		assert subtype in list(network['Subtype']), KeyError(f"{subtype} subtype not found in network layer of {name} tile")

		self.buildings = buildings[buildings['Subtype'] == subtype]
		self.parcels = parcels[parcels['Subtype'] == subtype]
		self.network = network[network['Subtype'] == subtype]

	def move_all_layers(self, reference):

		self.buildings = Shape(self.buildings).move(reference)
		self.parcels = Shape(self.parcels).move(reference)
		self.network = Shape(self.network).move(reference)
		return
