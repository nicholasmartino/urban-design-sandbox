from morphology.ShapeTools import Analyst
import geopandas as gpd


grid_gdf = gpd.read_file('/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Data/gmm_grids/fishnet_CoV_gmm_r4_built.shp')
directory = '/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Sandbox/shp/Tiers/GrossNewPopulation/Shp'
exp_cols = ['ToC2', 'CoN2', 'EcoCity2']

for scenario in ['Scenario1-GrossNEWpop.shp', 'Scenario2-GrossNEWpop.shp', 'Scenario3-GrossNEWpop.shp']:
	gdf = gpd.read_file(f'{directory}/{scenario}')
	out_cols = ['geometry']
	for experiment in exp_cols:
		for unique in gdf[experiment].unique():
			gdf.loc[gdf[experiment] == unique, f'{experiment}_{unique}'] = 1
			out_cols.append(f'{experiment}_{unique}')
	grid_gdf = Analyst(grid_gdf, gdf.loc[:, out_cols]).spatial_join(operations=['sum'])
grid_gdf.to_file(f'{directory}/fishned_CoV_gmm_joined.geojson', driver='GeoJSON')
