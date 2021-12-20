import pandas as pd
import geopandas as gpd

DIRECTORY = '/Volumes/SALA/Research/eLabs/50_projects/22_ESRI_Dashboard/Data'
LAYERS = {
	'Prince George.gdb': {
		'E0': 's_200326_jk_udes_E0_2020_BLDGS',
		'E1': 's_200316_jk_E7_2050_PPOLICY_BLDGS_tech_shell',
		'E2': 's_200316_jk_E5_2050_NC_FOCUS_BLDGS_tech_shell',
		'E3': 's_200316_jk_E6_2050_COR_FOCUS_BLDGS_tech_shell'
	},
	'Vancouver.gdb': {
		'E0': 'E0_2020_bldgs',
		# 'E1': 'E1_2030_D_bldgs',
		# 'E2': 'E2_2030_COR_A_bldgs',
		'E3': 'E3_2040_D_bldgs',
		# 'E4': 'E4_2040_COR_A_bldgs',
		'E5': 'E5_2040_COR_B_bldgs',
		'E6': 'E6_2040_TOD_B_bldgs',
		# 'E7': 'E7_2050_COR_B_bldgs',
		# 'E8': 'E8_2050_TOD_B_bldgs'
	},
	'Victoria.gdb': {
		'E0': 's_jk_200717_HQ_E0_2020_BLDGS_1',
		'E1': 's_jk_200808_HQ_E1_2040_BLDGS_1',
		'E2': 's_jk_200816_HQ_E2_2040_BLDGS_1',
		'E3': 's_jk_200808_HQ_E3_2040_BLDGS_1'
	}
}
COLUMNS = {
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
	'height': ['Height', 'Building_Height', 'height', 'BuidlingHeight', 'building_h'],
	'laneway': ['laneway'],
	'age': ['age', 'AGE'],
	'em_ht_t': ['em_ht_t'],
	'em_cl_t': ['em_cl_t'],
	'em_eqp_t': ['em_eqp_t'],
	'em_dhw_t': ['em_dhw_t'],
	'em_lgt_t': ['em_lgt_t'],
	'ann_em_t': ['ann_em_t']
}
NETWORK_LAYERS = {
	'Prince George.gdb': {
		'E0': {
			'Roads': 's_200107_cc_E0_2020_ROADS',
			'Transit': 's_200214_cc_E0_2020_TRANSIT'
		},
		'E1': {
			'Roads': 's_200214_cc_AT_2050_ROADS',
			'Transit': 's_200225_cc_AT_2050_TRANSIT',
		},
		'E2': {
			'Roads': 's_200214_cc_AT_2050_ROADS',
			'Transit': 's_200225_cc_AT_2050_TRANSIT',
		},
		'E3': {
			'Roads': 's_200214_cc_AT_2050_ROADS',
			'Transit': 's_200225_cc_AT_2050_TRANSIT',
		}
	},
	'Vancouver.gdb': {
		'E0': {
			'Roads': '',
			'Transit': ''
		},
		'E3': {
			'Roads': '',
			'Transit': ''
		},
		'E5': {
			'Roads': '',
			'Transit': ''
		},
		'E6': {
			'Roads': '',
			'Transit': ''
		},
	},
	'Victoria.gdb': {
		'E0': {
			'Roads': 'nm_200612_HQ_E0_2020_STRTS',
			'Transit': 'nm_200612_HQ_E0_2020_INTRS'
		},
		'E1': {
			'Roads': 'nm_200612_HQ_E0_2020_STRTS',
			'Transit': 'nm_200612_HQ_E0_2020_INTRS'
		},
		'E2': {
			'Roads': 'nm_200612_HQ_E0_2020_STRTS',
			'Transit': 'nm_200612_HQ_E0_2020_INTRS'
		},
		'E3': {
			'Roads': 'nm_200612_HQ_E0_2020_STRTS',
			'Transit': 'nm_200612_HQ_E0_2020_INTRS'
		}
	}
}


def get_d2use(g, land_use):
	return [geom.distance(g[g['landuse'].isin(land_use)].unary_union) for geom in g['geometry']]


if __name__ == '__main__':
	gdfs = gpd.GeoDataFrame()

	for city, experiments in LAYERS.items():
		for experiment, layer in experiments.items():
			gdf = gpd.read_file(f'{DIRECTORY}/{city}', layer=layer)
			gdf['city'] = city
			gdf['experiment'] = experiment

			# Estimate distances to cycling lanes
			if NETWORK_LAYERS[city][experiment]['Roads'] != '':
				roads_gdf = gpd.read_file(f'{DIRECTORY}/{city}', layer=NETWORK_LAYERS[city][experiment]['Roads'])
				if 'Bikeways' in roads_gdf.columns:
					bike_uu = roads_gdf[roads_gdf['Bikeways'] == 1].unary_union
				elif ('cycle_ocp' in roads_gdf.columns) & (experiment == 'E0'):
					bike_uu = roads_gdf[roads_gdf['cycle_ocp'] == 'existing'].unary_union
				else:
					bike_uu = roads_gdf[roads_gdf['cycle_2040'] == 1].unary_union
				gdf['d2bk'] = [geom.distance(bike_uu) for geom in gdf['geometry']]
			else:
				print(f"Distance to cycling lanes not calculated for {experiment} of {city}")

			# Estimate distances to transit stops
			if NETWORK_LAYERS[city][experiment]['Transit'] != '':
				transit_gdf = gpd.read_file(f'{DIRECTORY}/{city}', layer=NETWORK_LAYERS[city][experiment]['Transit'])
				if ('bus_2020' in transit_gdf.columns) & (experiment == 'E0'):
					transit_uu = transit_gdf[transit_gdf['bus_2020'] == 1].unary_union
				elif 'freqt_2040' in transit_gdf.columns:
					transit_uu = transit_gdf[transit_gdf['freqt_2040'] == 1].unary_union
				else:
					transit_uu = transit_gdf.unary_union
				gdf['d2tr'] = [geom.distance(transit_uu) for geom in gdf['geometry']]
			else:
				print(f"Distance to transit stops not calculated for {experiment} of {city}")

			# Transform emissions kg data into tons
			for em_col in ['em_ht_t', 'em_cl_t', 'em_eqp_t', 'em_dhw_t', 'em_lgt_t', 'ann_em_t']:
				if em_col not in gdf.columns:
					kg_col = f"{em_col.split('_t')[0]}_kg"
					assert kg_col in gdf.columns, AssertionError(f"No {em_col} or {kg_col} data for {layer} layer in {city} database")
					gdf[em_col] = gdf[kg_col] / 1000

			# Check if columns exist
			for new, old_columns in COLUMNS.items():
				common_cols = set(old_columns).intersection(set(gdf.columns))
				assert len(common_cols) > 0, AssertionError(f"No {new} data for {layer} layer in {city} database")
				for old in old_columns:
					if old in gdf.columns:
						gdf[new] = gdf[old]

			# Calculate distance indicators
			gdf['d2os'] = get_d2use(gdf, ['OS'])
			gdf['d2cv'] = get_d2use(gdf, ['CV'])
			gdf['d2cm'] = get_d2use(gdf, ['CM', 'MX'])

			gdfs = pd.concat([gdfs, gdf])
	gdfs = gdfs.loc[:, list(COLUMNS.keys()) + ['d2os', 'd2cv', 'd2cm', 'city', 'experiment', 'geometry']]

	# Rename land use
	rename = {
		'RS_SF_D': 'SFD',
		'RS_MF_L': 'MFL',
		'RS_SF_A': 'SFA',
		'RS_MF_H': 'MFH'
	}
	gdfs = gdfs.replace(rename)

	# Check for null values
	for col in gdfs.columns:
		nulls = gdfs[col].isna().sum()
		if nulls > 0:
			print(f"{nulls} null values in {col} column")
	gdfs.to_file(f"{DIRECTORY}/buildings_20211209.shp", driver='ESRI Shapefile')
