import base64
import io
import os.path
import sys
import time
import traceback

import dash
import geopandas as gpd
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
import pydeck as pdk
from dash import html
from dash.dependencies import Output, Input, State
from dash_deck import DeckGL
from shapely.geometry import Polygon

from layout.layout import layout
from layout.layout import types
from models.Grid import Grid
from models.Sandbox import Scenario, Indicators
from models.TestUploads import TestGridUpload
from models.Tiles import TILES
from store import TYPES, GRID_GDF, TILE_GDF, GRID_FILE, STREETS

external_stylesheets = [
    'https://fonts.googleapis.com/css2?family=Roboto&display=swap',
    'https://codepen.io/chriddyp/pen/bWLwgP.css'
]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.layout = layout

application = app.server
app.config.suppress_callback_exceptions = True

random_seed = 0
transparent = 'rgba(255,255,254,0)'

COLORS = {
    "SFD": [255, 248, 165],
    "SFA": [255, 237, 33],
    "MFL": [164, 126, 0],
    "MFM": [176, 173, 0],
    "MFH": [92, 76, 0],
    "MX": [255, 127, 0],
    "IND": [128, 128, 128],
    "CM": [200, 90, 90],
    "CV": [100, 172, 190],
    "OS": [180, 210, 180],
}
color_discrete_map = {key: f'rgb{tuple(i)}' for key, i in COLORS.items()}
template = dict(layout=go.Layout(
    title_font=dict(family="Roboto Black", size=16, color='DarkGray'),
    font=dict(family="Roboto Light"),
    margin=dict(l=60, r=20, t=40, b=60),
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)',
    showlegend=False,
))
# Dash Leaflet parameters
MAP_ID = "map-id"
POLYLINE_ID = "polyline-id"
POLYGON_ID = "polygon-id"
dummy_pos = [0, 0]
dlatlon2 = 1e-6  # Controls tolerance of closing click
DIRECTORY = \
    "/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership " \
    "Engage/Sandbox/shp/Tiers/GrossNewPopulation/Shp/ "
all_tiles = TILE_GDF.copy()

# Define zoning layer colors

# Get open spaces from CoV open data
if sys.platform == 'win32':
    GPKG = 'G:/My Drive/Databases/Vancouver, British Columbia.gpkg'
else:
    GPKG = '/Volumes/Samsung_T5/Databases/Vancouver, British Columbia.gpkg'
PARKS = gpd.read_feather('data/feather/cov_parks.feather')
# REAL_TREES = gpd.read_feather('data/feather/nature_tree_canopy.feather')
land_use_gdf = gpd.read_feather('data/feather/mvan_landuse.feather')
diagonal_gdf = gpd.read_file('data/geojson/diagonal_tiles.geojson')
CLUSTERS = {
    2: 'cell_types8.png',
    3: 'cell_types3.png',
    4: 'cell_types7.png',
    5: 'cell_types2.png',
    8: 'cell_types6.png',
    9: 'cell_types.png',
    10: 'cell_types4.png',
    11: 'cell_types5.png',
    12: 'cell_types.png',
}


def gdf_from_memory(memory):
    geom_types = [i.__class__.__name__ for i in memory['object']['geometry']]
    if 'Polygon' not in geom_types:
        try:
            memory['object']['geometry'] = [Polygon(coord) for coord in memory['object']['geometry']['coordinates']]
        except:
            memory['object']['geometry'] = [Polygon(coord) for coord in memory['object']['geometry']['coordinates'][0]]
    gdf = gpd.GeoDataFrame(memory['object'], geometry='geometry', crs=4326)
    if 'id' not in gdf.columns:
        gdf['id'] = gdf.reset_index(drop=True).index
    return gdf


def pick_geometry(deck, gdf, cell_id):
    new_layer = pdk.Layer(
        id=cell_id,
        type="GeoJsonLayer",
        data=gdf.to_crs(4326),
        extruded=False,
        opacity=0.5,
        getFillColor=[0, 0, 255],
        get_line_color=[255, 255, 255],
        auto_highlight=True,
        pickable=True,
    )
    if cell_id not in [lay.id for lay in deck.layers]:
        deck.layers.append(new_layer)
    else:
        deck.layers = [lay for lay in deck.layers if lay.id != cell_id]
    return deck


def create_bld_layer(l_use, colors, buildings):
    return pdk.Layer(
        id=l_use,
        type="GeoJsonLayer",
        extruded=True,
        getElevation="height",
        opacity=0.8,
        data=buildings.loc[buildings['LANDUSE'] == l_use, ['height', 'geometry']].to_crs(4326),
        getFillColor=colors[l_use],
    )


def create_pcl_layers(gdf):
    parcels_pdk = pdk.Layer(
        id=f"parcels",
        type="GeoJsonLayer",
        opacity=0.2,
        data=gdf.loc[gdf['LANDUSE'] != 'OS', ['geometry']].to_crs(4326),
        getFillColor=[180, 180, 180],
    )
    open_pdk = pdk.Layer(
        id=f"open_spaces",
        type="GeoJsonLayer",
        opacity=0.5,
        data=gdf.loc[gdf['LANDUSE'] == 'OS', ['geometry']].to_crs(4326),
        getFillColor=COLORS['OS'],
    )

    return parcels_pdk, open_pdk


def build_scenario(tiles, prefix):
    prcls = tiles[tiles['Type'] == 'prcls']
    bldgs = tiles[tiles['Type'] == 'bldgs']
    strts = tiles[tiles['Type'] == 'strts']
    blcks = tiles[tiles['Type'] == 'block']
    trees = tiles[tiles['Type'] == 'trees']

    scn = Scenario(
        parcels=prcls, buildings=bldgs, trees=trees, real_parks=PARKS, real_trees=None, name=prefix
    )
    scn.parcels = scn.extract_parks()
    # sb_trees = scn.extract_trees(directory='data/feather')
    # sb_trees.to_feather(f"data/feather/{prefix}trees.feather")

    ind = Indicators(parcels=scn.parcels, buildings=bldgs, streets=strts, blocks=blcks)
    ind.test_indicators()

    print(f"Total population: {ind.get_total_population()}")
    print(f"Total area: {ind.get_total_area()}")
    return ind


def join_high_st(gdf):
    STREETS['geometry'] = STREETS.buffer(5)
    gdf['id'] = gdf.index
    gdf.loc[gpd.overlay(gdf, STREETS[STREETS['Category'] == 'Arterial'])['id'], 'Arterial'] = 1
    gdf.loc[gdf['Type'].isin(['Mid_High_Street', 'Moderate_Density', 'Dense_Nodal']), 'High St Type'] = 1
    gdf.loc[(gdf['Arterial'] == 1) & (gdf['High St Type'] == 1), 'High St'] = 1
    return gdf


def parse_geojson(encoded):
    content_type, content_string = encoded.split(',')
    decoded = base64.b64decode(content_string)
    return gpd.read_file(io.BytesIO(decoded)).to_crs(26910)


# Define initial view state and create deck object
view_state = pdk.ViewState(latitude=49.254, longitude=-123.13, zoom=15, max_zoom=20, pitch=60, bearing=0)
r = pdk.Deck(
    layers=[],
    initial_view_state=view_state,
    map_style=pdk.map_styles.LIGHT_NO_LABELS,
)
mapbox_key = "pk.eyJ1IjoibmljaG9sYXNtYXJ0aW5vIiwiYSI6ImNrMjVhOGphOTAzZGUzbG8wNHJhdTZrMmYifQ.98uDMnGIvn1zrw4ZWUO35g"

tooltip = {
    "html": "F.S.R. {Max FSR} <br> Height {Max Height}"
}

html.Button(id='rotate', value='Rotate'),


# STORES
def assign_callback(app_, event):
    @app_.callback(
        Output(component_id=f"memory", component_property="data"),
        [Input(component_id="deck", component_property=event)],
    )
    def dump_json(data):
        return data


assign_callback(app, "clickInfo")


@app.callback(
    Output(component_id="cell_type_icon", component_property="figure"),
    Output(component_id="selected_type", component_property="children"),
    [Input(component_id="memory", component_property="data")]
)
def update_image(memory):
    print(update_image)

    file_name = CLUSTERS[memory['object']['clus_gmm']]
    type_name = TYPES[memory['object']['clus_gmm']]

    fig = go.Figure()

    fig.add_layout_image(
        dict(
            source=f"https://raw.githubusercontent.com/nicholasmartino/urban-design-sandbox/master/images/{file_name}",
            xref="x", x=0, sizex=8,
            yref="y", y=5, sizey=8,
            opacity=0.8,
            layer="below"
        )
    )
    fig.update_layout(template=template, margin=dict(l=0, r=0, t=0, b=0))
    fig.update_xaxes(visible=False, showticklabels=False)
    fig.update_yaxes(visible=False, showticklabels=False)
    return fig, type_name


# LAYERS + DESIGN
@app.callback(
    Output(component_id="deck_div", component_property="children"),
    Output(component_id="test-error", component_property="displayed"),
    Output(component_id="test-error", component_property="message"),
    Output(component_id="area_by_lu", component_property="figure"),
    Output(component_id="fsr_hist", component_property="figure"),
    Output(component_id="dwelling_mix", component_property="figure"),
    Output(component_id="proximity", component_property="figure"),
    Output(component_id="total_units", component_property="children"),
    Output(component_id="total_population", component_property="children"),
    Output(component_id="fsr", component_property="children"),
    Output(component_id="max_height", component_property="children"),
    Output(component_id="loading-output-1", component_property="children"),
    [
        # Input(component_id="type", component_property="value"),
        # Input(component_id='rotate', component_property='n_clicks'),
        # Input(component_id='flip_h', component_property='n_clicks'),
        # Input(component_id='flip_v', component_property='n_clicks'),
        # Input(component_id='select', component_property='n_clicks'),
        Input(component_id='update_sandbox', component_property='data'),
        Input(component_id='memory', component_property='data'),
        Input(component_id='upload', component_property='contents'),
        State(component_id='upload', component_property='filename'),
    ])  # , State('input-on-submit', 'value'))
def main_callback(update_sandbox, memory, uploaded, file_name):
    stt = time.time()
    types_named = {t: i for i, t in enumerate(types)}

    if uploaded is not None:
        grid_gdf = parse_geojson(uploaded)
        try:
            TestGridUpload(grid_gdf).test_all()
        except AssertionError:
            _, _, tb = sys.exc_info()
            return None, True, f"{_}", None, None, None, None, None, None, None, None, None
        prefix = file_name
        r.layers = []
    else:
        grid_gdf = GRID_GDF
        prefix = GRID_FILE.split('.')[0]

    all_layers_file = f"data/feather/{prefix}_all_layers.feather"
    buildings_file = f"data/feather/{prefix}buildings.feather"
    parcels_file = f"data/feather/{prefix}parcels.feather"

    grid_gdf['id'] = grid_gdf.index

    # Get callback context
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = 'No clicks yet'
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # Create cells deck layer
    if 'cells' not in [l.id for l in r.layers]:
        cells = pdk.Layer(
            id=f"cells",
            type="GeoJsonLayer",
            data=grid_gdf.to_crs(4326),
            stroked=True,
            wireframe=True,
            opacity=0.2,
            getFillColor=[180, 180, 180, 255],
            get_line_color=[0, 0, 0],
            get_line_width=4,
            auto_highlight=True,
            pickable=True,
        )
        r.layers.append(cells)

    # Append selected parcel if exists
    if button_id == 'memory':
        cell_gdf = gdf_from_memory(memory)
        cell_id = f"{cell_gdf.loc[0, 'id']}"
        pick_geometry(r, cell_gdf, cell_id)

    if memory is None:
        # Run grids
        grid_gdf['id'] = grid_gdf.index
        grid_gdf['Type'] = grid_gdf['clus_gmm'].replace(TYPES)
        grid_gdf.loc[
            grid_gdf['Type'].isin(['Mid_High_Street', 'Moderate_Density', 'Dense_Nodal', 'Dense_TOD',
                                   'Dense_Industrial']), 'High St Type'] = 1
        if 'High St' not in grid_gdf.columns:
            grid_gdf = join_high_st(grid_gdf)

        # Read/export all layers
        if (os.path.exists(all_layers_file)) \
                and (os.path.exists(buildings_file)) \
                and (os.path.exists(parcels_file)) \
                and (not update_sandbox):
            tiles = gpd.read_feather(all_layers_file)
        else:
            grid = Grid(grid_gdf, TILES, prefix=f"{prefix}_", land_use=land_use_gdf, diagonal_gdf=diagonal_gdf)
            grid.test_assign_subtypes()
            tiles = grid.test_place_tiles()
            tiles.loc[:, [col for col in tiles.columns if col not in ['laneway']]].to_feather(all_layers_file)

        prefix = f"{prefix.split('.')[0]}_"

        p_exists = os.path.exists(f"data/feather/{prefix}parcels.feather")
        b_exists = os.path.exists(f"data/feather/{prefix}buildings.feather")
        # t_exists = os.path.exists(f"data/feather/{prefix}trees.feather")

        if p_exists and b_exists:
            parcels = gpd.read_feather(f'data/feather/{prefix}parcels.feather')
            buildings = gpd.read_feather(f'data/feather/{prefix}buildings.feather')
            # trees = gpd.read_feather(f'data/feather/{prefix}trees.feather')
            ind = Indicators(parcels=parcels, buildings=buildings)

        else:
            print(f"\n {prefix}")
            tiles = tiles.reset_index(drop=True)
            ind = build_scenario(tiles, prefix)
            parcels = ind.parcels.copy()
            buildings = ind.buildings.copy()
            ind.parcels.to_feather(f"data/feather/{prefix}parcels.feather")
            ind.buildings.loc[:, [c for c in ind.buildings.columns if c not in ['comm_units']]]. \
                to_feather(f"data/feather/{prefix}buildings.feather")
            tiles.loc[:, [col for col in tiles.columns if col not in ['laneway']]].to_feather(all_layers_file)

        # Create buildings and parcels layers
        parcels_pdk, open_pdk = create_pcl_layers(parcels)
        r.layers = r.layers + [parcels_pdk, open_pdk]
        for use in COLORS.keys():
            r.layers.append(create_bld_layer(use, COLORS, buildings))

        # Land use
        area_by_lu = px.bar(
            ind.get_floor_area_by_land_use(), x='Land Use', y='Floor Area (m²)', title='FLOOR AREA BY LAND USE',
            template=template, color='Land Use', color_discrete_map=color_discrete_map
        )

        # Floor space ratio
        fsr = ind.get_parcel_far().dropna(subset=['far'])
        res_uses = ['SFD', 'SFA', 'MFM', 'MFL', 'MFH']
        fsr.loc[fsr['Land Use'].isin(res_uses), 'lu_type'] = 'Residential'
        fsr.loc[~fsr['Land Use'].isin(res_uses), 'lu_type'] = 'Non-Residential'
        fsr.loc[fsr['Land Use'] == 'MX', 'lu_type'] = 'Mixed'
        fsr_hist = ff.create_distplot(
            [list(fsr.loc[fsr['lu_type'] == u, 'Floor Area Ratio']) for u in fsr['lu_type'].unique()],
            group_labels=fsr['lu_type'].unique(), colors=[
                f"rgb{tuple(COLORS['CM'])}", f"rgb{tuple(COLORS['MX'])}", f"rgb{tuple(COLORS['MFM'])}"
            ]
        )
        fsr_hist.update_layout(
            template['layout'],
            title={'text': 'FLOOR / PARCEL AREA RATIO (FAR)'},
            xaxis_title='FAR',
            yaxis_title='Count (%)',
        )

        # Dwelling mix
        dwelling_mix = px.bar(
            ind.get_dwelling_mix(), x='Floor Area', y='Dwelling Type', template=template,
            color='Dwelling Type', color_discrete_map=color_discrete_map, title='DWELLING MIX'
        )

        # Proximity
        prx_df = ind.get_proximity()
        proximity = px.scatter(
            prx_df, y="Land Use", x="Proximity (m)", color="Land Use", title='PROXIMITY',
            color_discrete_map=color_discrete_map, template=template
        )
        proximity.update_traces(marker=dict(opacity=0.1))
        for use in ['OS', 'CV', 'MX', 'CM']:
            mean = prx_df.loc[prx_df['land_use'] == use, 'proximity'].mean()
            proximity.add_shape(y0=use, y1=use, x0=mean, x1=mean + 10, fillcolor='black')
    else:
        prefix = GRID_FILE.split('.')[0]
        prefix = f"{prefix.split('.')[0]}_"
        parcels = gpd.read_feather(f'data/feather/{prefix}parcels.feather')
        buildings = gpd.read_feather(f'data/feather/{prefix}buildings.feather')

        cells = gdf_from_memory(memory).to_crs(26910)
        cells['Type'] = cells['clus_gmm'].replace(TYPES)

        in_cell_bld = gpd.overlay(buildings, cells.loc[:, ['geometry']])
        in_cell_pcl = gpd.overlay(parcels, cells.loc[:, ['geometry']])

        ind = Indicators(parcels=in_cell_pcl, buildings=in_cell_bld)
        area_by_lu = px.bar(
            ind.get_area_by_land_use(), x='Land Use', y='Area (m²)', template=template,
            color='Land Use', color_discrete_map=color_discrete_map
        )

    # Create deck-gl object
    dgl = DeckGL(
        id="deck",
        data=r.to_json(),
        mapboxKey=mapbox_key,
        enableEvents=['click'],
        style={'width': '100%', 'float': 'left', 'display': 'inline-block'},
    )

    total_units = int(ind.get_residential_units()['res_units'].sum())
    total_population = int(ind.get_resident_count()['res_count'].sum())
    fsr = round(sum(ind.get_floor_area_by_land_use()['Floor Area (m²)']) / sum(
        ind.parcels[ind.parcels['LANDUSE'] != 'OS'].area), 2)
    max_height = f"Max height: {max(ind.buildings['height'])} m ({int(max(ind.buildings['height']) / 3)} stories)"
    print(f"Callback: {round((time.time() - stt), 3)} seconds with {[l.id for l in r.layers]} layers")
    return \
        dgl, None, None, area_by_lu, fsr_hist, dwelling_mix, proximity, f"{total_units} units", \
        f"{total_population} people", f"Mean FSR: {fsr}", max_height, None


# Download contents
@app.callback(
    Output("download", "data"),
    Input("deck_div", "children"),
    Input("btn-download", "n_clicks"),
    prevent_initial_call=True,
)
def download_layers(deck_div, n_clicks):
    if n_clicks is not None:
        out_gdf = gpd.GeoDataFrame()
        for file in os.listdir('data/feather'):
            if 'all_layers' in file:
                out_gdf = pd.concat([out_gdf, gpd.read_feather(f'data/feather/{file}')])
        return dict(content=out_gdf.to_json(), filename="all_layers.geojson")


if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8080)  # Run in the elementslab server
