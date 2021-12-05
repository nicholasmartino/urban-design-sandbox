import geopandas as gpd

grid_gdf = gpd.read_file('/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Data/gmm_grids/fishnet_CoV_gmm_r4_built.shp')
directory = '/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Sandbox/shp/Tiers/GrossNewPopulation/Shp'
exp_cols = ['ToC2', 'CoN2', 'EcoCity2']

grid_gdf['grid_id'] = grid_gdf.index
for scenario in ['Scenario1-GrossNEWpop.shp', 'Scenario2-GrossNEWpop.shp', 'Scenario3-GrossNEWpop.shp']:
	gdf = gpd.read_file(f'{directory}/{scenario}')
	gdf['id'] = gdf.index
	out_cols = ['geometry']
	for experiment in exp_cols:
		joined_grid = gpd.sjoin(grid_gdf, gdf.loc[:, ['id', experiment, 'geometry']])
		joined_grid = joined_grid.dropna(subset=[experiment])
		grouped = joined_grid.groupby('grid_id', as_index=False).first()
		grouped.index = grouped['grid_id']
		grid_gdf.loc[grouped['grid_id'], f'{scenario}_{experiment}'] = grouped[experiment]

grid_gdf.to_file(f'{directory}/fishned_CoV_gmm_joined.geojson', driver='GeoJSON')
