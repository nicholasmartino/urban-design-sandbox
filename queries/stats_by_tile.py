import os

import geopandas as gpd
import pandas as pd
os.chdir('/Volumes/GoogleDrive/My Drive/Python/urban-design-sandbox')

from Inputs import SERVER_DIR
from Sandbox import Indicators

TILES = gpd.read_feather('data/feather/all_tiles.feather')

summary = pd.DataFrame()
for cell_type in TILES['cell_type'].unique():
	for subtype in TILES['Subtype'].unique():
		st_tiles = TILES[(TILES['cell_type'] == cell_type) & (TILES['Subtype'] == subtype)]
		if len(st_tiles) > 0:
			parcels = st_tiles[st_tiles['Type'] == 'prcls']
			buildings = st_tiles[st_tiles['Type'] == 'bldgs']

			ind = Indicators(parcels=parcels, buildings=buildings)
			ind.buildings = ind.get_buildings_floor_area()
			ind.buildings = ind.get_residential_units()
			ind.buildings = ind.get_commercial_units()
			ind.buildings = ind.get_resident_count()
			ind.buildings = ind.get_height_from_storeis()

			i = len(summary)
			summary.loc[i, 'Type'] = cell_type
			summary.loc[i, 'Subtype'] = subtype
			summary.loc[i, 'Residents'] = ind.buildings['res_count'].sum()
			summary.loc[i, 'FAR'] = round(sum(ind.buildings['floor_area'])/sum(ind.parcels.area), 2)
			summary.loc[i, 'Max Height (m)'] = ind.buildings['height'].max()

			area_by_lu = ind.get_area_by_land_use()
			for land_use in area_by_lu['Land Use'].unique():
				summary.loc[i, f'{land_use} (m²)'] = round(sum(
					area_by_lu.loc[(area_by_lu['Land Use'] == land_use), 'Area (m²)']), 2)

summary = summary.fillna(0)
summary['Residents'] = summary['Residents'].astype(int)
summary['Max Height (m)'] = summary['Max Height (m)'].astype(int)
summary = summary.sort_values('FAR', ascending=False)
summary.to_csv(f'{SERVER_DIR}/Sandbox/Stats by Tiles.csv')
print(summary)
print(f'{SERVER_DIR}/Sandbox/Stats by Tiles.csv')
