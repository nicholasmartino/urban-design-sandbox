import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd

from variables import TYPES, COLORS

broadway_dir = f"S:/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Sandbox/shp/Broadway/Broadway Model Layers"

baseline_gdf = gpd.read_file(f"{broadway_dir}/Baseline_Fishnet_BroadwayPlan_ELABTypes_v3.shp")
plan_gdf = gpd.read_file(f"{broadway_dir}/Fishnet-BroadwayPlan-ELABTypes.shp")
plan_gdf = pd.merge(plan_gdf, baseline_gdf.loc[:, ['ID', 'cell_type', 'Subtype']], on="ID")


###
# Types
###

plan_gdf['clus_gmm'] = plan_gdf['cell_type'].replace(dict(map(reversed, TYPES.items())))
plan_gdf = plan_gdf.dropna(subset=['clus_gmm'])
plan_gdf.to_file('../data/geojson/broadway_baseline.geojson', driver='GeoJSON')

plan_gdf.loc[plan_gdf['ZoneType'] == 'Station Area', 'cell_type'] = "Dense_TOD"
plan_gdf.loc[plan_gdf['ZoneType'] == 'Shoulder Area', 'cell_type'] = "Moderate_Density"


plan_gdf.loc[
    (plan_gdf['ZoneType'] == 'Village'),  # Filter
    'cell_type'  # Column to update
] = "Mid_High_Street"  # Value to update

plan_gdf.loc[
    plan_gdf['ZoneType'] == 'Residential RM/FM',
    ('cell_type', 'Subtype')
] = ("Treed_Large_Home", "Multi-Family")

plan_gdf.loc[
    (plan_gdf['ZoneType'] == "Residential RT/RS Zone") &
    (plan_gdf['cell_type'] == "Van Detached"),
    'cell_type'
] = "Typical_Van_SF"

plan_gdf.loc[
    (plan_gdf['ZoneType'] == "Industrial/Employment Area") &
    (plan_gdf['cell_type'] == "Dense_Nodal"),
    "cell_type"
] = "Dense_Industrial"

plan_gdf.loc[
    (plan_gdf['ZoneType'] == "Industrial/Employment Area") &
    (~plan_gdf['cell_type'].isin(["Dense_Nodal", "Dense_Industrial"])),
    "cell_type"
] = "Moderate_Industrial"


plan_gdf = plan_gdf[~plan_gdf['cell_type'].isna()]
plan_gdf['clus_gmm'] = plan_gdf['cell_type'].replace(dict(map(reversed, TYPES.items())))
plan_gdf = plan_gdf.dropna(subset=['clus_gmm'])
plan_gdf.to_file('../data/geojson/broadway_future.geojson', driver='GeoJSON')
