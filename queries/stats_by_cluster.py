import geopandas as gpd

MV = gpd.read_file('/Users/nicholasmartino/Desktop/Landuse2016/Landuse2016.shp')
CL = gpd.read_file('../data/clustering_summary.geojson')

uses = list([i for i in MV['elab_lu'].unique() if i is not None])
for cluster in CL['clus_gmm'].unique():
	cluster_gdf = CL[CL['clus_gmm'] == cluster]
	for use in uses:
		overlay = gpd.overlay(cluster_gdf, MV.loc[MV['elab_lu'] == use, ['geometry']])
		CL.loc[CL['clus_gmm'] == cluster, f'area_{use}'] = overlay.area.sum()

for cluster in CL['clus_gmm'].unique():
	cluster_gdf = CL[CL['clus_gmm'] == cluster]
	for use in uses:
		CL.loc[CL['clus_gmm'] == cluster, f'share_{use}'] = round(cluster_gdf[f'area_{use}'].sum()/sum(cluster_gdf.loc[:, [f'area_{i}' for i in uses]].sum()), 2)

CL.to_crs(4326).to_file('../data/clustering_summary_lu.geojson', driver='GeoJSON')
