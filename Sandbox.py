import pandas as pd
import geopandas as gpd
from morphology.ShapeTools import Analyst


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

	def get_residential_units(self, unit_area=75):
		gdf = self.buildings.copy()
		gdf.loc[gdf['LANDUSE'].isin(['CM', 'CV', 'IND', 'OS']), 'res_units'] = 0
		gdf.loc[gdf['LANDUSE'].isin(['MFL', 'MFM', 'MFH']), 'res_units'] = (gdf['floor_area']/unit_area).astype(int) * gdf['maxstories']
		gdf.loc[gdf['LANDUSE'] == 'MX', 'res_units'] = (gdf['floor_area']/unit_area).astype(int) * (gdf['maxstories'] - 1)
		gdf['res_units'] = gdf['res_units'].fillna(0)
		return gdf

	def get_commercial_units(self, unit_area=100):
		gdf = self.buildings.copy()
		gdf.loc[gdf['LANDUSE'] == 'CM', 'comm_units'] = (gdf['floor_area']/unit_area).astype(int) * gdf['maxstories']
		gdf.loc[gdf['LANDUSE'] == 'MX', 'comm_units'] = (gdf['floor_area']/unit_area).astype(int)
		return gdf

	def get_resident_count(self):
		gdf = self.buildings.copy()
		gdf['res_count'] = gdf['res_units'].astype(int) * 1.82
		return gdf

	def remove_buildings_from_open_spaces(self):
		gdf = self.buildings.copy()
		gdf['geometry'] = gdf.centroid
		os_bld = gpd.overlay(gdf, self.parcels[self.parcels['LANDUSE'] == 'OS'])
		return self.buildings[~self.buildings['bid'].isin(os_bld['bid'])]

	def get_area_by_land_use(self):
		gdf = self.parcels.copy()
		return pd.DataFrame(gdf.groupby('LANDUSE').sum()['area'])

	def get_population_stats(self):
		gdf = self.buildings.copy()
		return gdf['res_count'].sum()

	def test_indicators(self):
		bld_cols = ['LANDUSE', 'maxstories', 'laneway']
		for col in bld_cols:
			assert col in self.buildings.columns, KeyError(f'{col} column not found in buildings layer')
		self.buildings['maxstories'] = pd.to_numeric(self.buildings['maxstories'])
		assert 'LANDUSE' in self.parcels.columns, KeyError('LANDUSE column not found in parcels layer')

		self.buildings = self.get_buildings_floor_area()
		self.parcels = self.get_parcel_far()
		self.buildings = self.get_residential_units()
		self.buildings = self.get_commercial_units()
		self.buildings = self.get_resident_count()
		self.buildings = self.remove_buildings_from_open_spaces()
		self.get_area_by_land_use()

		return


if __name__ == '__main__':
	ind = Indicators(parcels=gpd.read_file('data/parcels.shp'), buildings=gpd.read_file('data/buildings.shp'))
	parks = gpd.read_file('https://opendata.vancouver.ca/explore/dataset/parks-polygon-representation/download/?format=geojson&timezone=America/Los_Angeles&lang=en&').to_crs(26910)
	pcl = ind.parcels.copy()
	pcl['geometry'] = pcl.centroid
	inters_parcels = gpd.overlay(pcl.loc[:, ['pid', 'geometry']], parks)
	ind.parcels.loc[ind.parcels['pid'].isin(inters_parcels['pid']), 'LANDUSE'] = 'OS'
	ind.test_indicators()
	ind.parcels.to_file('data/parcels_indicator.shp')
	ind.buildings.to_file('data/buildings_indicator.shp')
	ind.get_area_by_land_use().to_csv('data/land_use.csv')
