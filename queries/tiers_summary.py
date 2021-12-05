import geopandas as gpd
import matplotlib.pyplot as plt

from Inputs import TYPES

directory = '/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Sandbox/shp/Tiers/GrossNewPopulation/Shp'

# Plot the proportion of tiers in baseline elementslab and future scenarios from CoV
grid_gdf = gpd.read_file(f'{directory}/fishned_CoV_gmm_joined.geojson', driver='GeoJSON')
grid_gdf['clus_gmm_names'] = grid_gdf['clus_gmm'].replace(TYPES)
to_plot = ['clus_gmm_names', 'Scenario1-GrossNEWpop.shp_ToC2', 'Scenario1-GrossNEWpop.shp_CoN2', 'Scenario1-GrossNEWpop.shp_EcoCity2']

fig, ax = plt.subplots(1, len(to_plot), figsize=(15, 8))
for i, col in enumerate(to_plot):
	grid_gdf[col].value_counts().plot(kind='bar', ax=ax[i])
	ax[i].set_title(col)
plt.gcf().subplots_adjust(bottom=0.25)
plt.savefig('tiers_summary.png')
