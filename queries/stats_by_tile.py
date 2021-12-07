import os

import geopandas as gpd
import pandas as pd

os.chdir('/Volumes/GoogleDrive/My Drive/Python/urban-design-sandbox')
from Sandbox import Indicators

TILES = gpd.read_feather('data/feather/all_tiles.feather')

summary = pd.DataFrame()
for cell_type in TILES['cell_type'].unique():
	st_tiles = TILES[TILES['cell_type'] == cell_type]
	parcels = st_tiles[st_tiles['Type'] == 'prcls']
	buildings = st_tiles[st_tiles['Type'] == 'bldgs']

	ind = Indicators(parcels=parcels, buildings=buildings)
	ind.buildings = ind.get_buildings_floor_area()
	ind.buildings = ind.get_residential_units()
	ind.buildings = ind.get_commercial_units()
	ind.buildings = ind.get_resident_count()
	ind.buildings = ind.get_height_from_storeis()

	summary.loc[cell_type, 'Residents'] = ind.buildings['res_count'].sum()/len(ind.parcels['Subtype'].unique())
	summary.loc[cell_type, 'FAR'] = round(sum(ind.buildings['floor_area'])/sum(ind.parcels.area), 2)
	summary.loc[cell_type, 'Max Height (m)'] = ind.buildings['height'].max()

	area_by_lu = ind.get_area_by_land_use()
	for land_use in area_by_lu['Land Use'].unique():
		summary.loc[cell_type, f'{land_use} (m²)'] = round(sum(
			area_by_lu.loc[(area_by_lu['Land Use'] == land_use), 'Area (m²)'])/len(ind.parcels['Subtype'].unique()), 2)

summary = summary.fillna(0)
summary['Residents'] = summary['Residents'].astype(int)
summary['Max Height (m)'] = summary['Max Height (m)'].astype(int)
summary = summary.sort_values('FAR', ascending=False)
print(summary)
summary.to_csv('/Users/nicholasmartino/Desktop/Stats by Tiles.csv')
