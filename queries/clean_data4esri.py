import pandas as pd
import geopandas as gpd

directory = '/Volumes/team/SALA/Research/eLabs/50_projects/22_ESRI_Dashboard/Data'
layers = {
	'Prince George.gdb': {
		'E0': 's_200326_jk_udes_E0_2020_BLDGS',
		'E1': 's_200316_jk_E7_2050_PPOLICY_BLDGS_tech_shell',
		'E2': 's_200316_jk_E5_2050_NC_FOCUS_BLDGS_tech_shell',
		'E3': 's_200316_jk_E6_2050_COR_FOCUS_BLDGS_tech_shell'
	},
	'Vancouver.gdb': {
		'E0': 'E0_2020_bldgs',
		'E1': 'E1_2030_D_bldgs',
		'E2': 'E2_2030_COR_A_bldgs',
		'E3': 'E3_2040_D_bldgs',
		'E4': 'E4_2040_COR_A_bldgs',
		'E5': 'E5_2040_COR_B_bldgs',
		'E6': 'E6_2040_TOD_B_bldgs',
		'E7': 'E7_2050_COR_B_bldgs',
		'E8': 'E8_2050_TOD_B_bldgs'
	},
	'Victoria.gdb': {
		'E0': 's_jk_200717_HQ_E0_2020_BLDGS_1',
		'E1': 's_jk_200808_HQ_E1_2040_BLDGS_1',
		'E2': 's_jk_200816_HQ_E2_2040_BLDGS_1',
		'E3': 's_jk_200808_HQ_E3_2040_BLDGS_1'
	}
}
columns = {
	'landuse': ['landuse', 'LANDUSE', 'Landuse', 'landuse_2'],
	'total_kwh': ['total_kWh'],
	'total_gj': ['total_GJ'],
	'eui_kwh': ['EUI_kWh'],
	'euc_kwh': ['EUC_kWh'],
	'heating_gj': ['heating_GJ_annual', 'heating_GJ'],
	'cooling_gj': ['cooling_GJ_annual', 'cooling_GJ'],
	'light_gj': ['lighting_GJ_annual', 'lighting_G', 'lighting_GJ'],
	'equip_gj': ['equipment_', 'equip_GJ', 'equipment_GJ_annual'],
	'dhw_gj': ['dhw_GJ_ann', 'dhw_GJ_annual'],
	'res_count': ['res_count'],
	'res_units': ['res_units', 'n_res_unit'],
	'floor_area': ['floor_area', 'TFA'],
	'mx_stories': ['maxstories'],
	'height': ['Height', 'Building_Height', 'height', 'BuidlingHeight', 'building_h']
}

if __name__ == '__main__':
	gdfs = gpd.GeoDataFrame()
	for city, experiments in layers.items():
		for experiment, layer in experiments.items():
			gdf = gpd.read_file(f'{directory}/{city}', layer=layer)
			gdf['city'] = city
			gdf['experiment'] = experiment
			for new, old_columns in columns.items():
				common_cols = set(old_columns).intersection(set(gdf.columns))
				assert len(common_cols) > 0, AssertionError(f"No {new} data for {layer} layer in {city} database")
				for old in old_columns:
					if old in gdf.columns:
						gdf[new] = gdf[old]
			gdfs = pd.concat([gdfs, gdf])
	gdfs = gdfs.loc[:, list(columns.keys()) + ['geometry']]

	# Rename land use
	rename = {
		'RS_SF_D': 'SFD',
		'RS_MF_L': 'MFL',
		'RS_SF_A': 'SFA',
		'RS_MF_H': 'MFH'
	}
	gdfs = gdfs.rename(rename)

	# Check for null values
	for col in gdfs.columns:
		nulls = gdfs[col].isna().sum()
		if nulls > 0:
			print(f"{nulls} null values in {col} column")
