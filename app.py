import os.path
import time

import dash
import dash_core_components as dcc
import dash_html_components as html
import numpy as np
import plotly.express as px
import pydeck as pdk
from dash.dependencies import Output, Input
from dash_bootstrap_components import Col
from dash_deck import DeckGL
from shapely.affinity import translate
from shapely.geometry import Polygon

from Grid import Grid
from Inputs import *
from Inputs import TYPES
from Sandbox import Scenario, Indicators
from Tiles import TILES

alpha = 100
random_seed = 0
transparent = 'rgba(255,255,255,0)'

# Dash Leaflet parameters
MAP_ID = "map-id"
POLYLINE_ID = "polyline-id"
POLYGON_ID = "polygon-id"
dummy_pos = [0, 0]
dlatlon2 = 1e-6  # Controls tolerance of closing click
DIRECTORY = "/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership Engage/Sandbox/shp/Tiers/GrossNewPopulation/Shp/"
GRID_FILE = "fishned_CoV_gmm_joined.geojson"
all_tiles = gpd.read_file(
	'/Volumes/SALA/Research/eLabs/50_projects/20_City_o_Vancouver/SSHRC Partnership '
	'Engage/Sandbox/shp/elementslab/Version_3/all_tiles.shp')

# Define zoning layer colors
types = {
	'Typical Van SF': [255, 100, 100, alpha],
	'2': [255, 100, 255, alpha],
	'3': [100, 100, 255, alpha],
	'4': [255, 255, 100, alpha],
	'5': [255, 127, 100, alpha],
	'6': [255, 127, 100, alpha],
	'7': [255, 127, 100, alpha],
	'8': [255, 127, 100, alpha],
	'9': [255, 127, 100, alpha],
	'10': [255, 127, 100, alpha],
	'11': [255, 127, 100, alpha],
}
grid_gdf = gpd.read_feather('data/feather/broadway_plan_baseline.feather')
# Get open spaces from CoV open data
if sys.platform == 'win32':
	GPKG = 'G:/My Drive/Databases/Vancouver, British Columbia.gpkg'
else:
	GPKG = '/Volumes/Samsung_T5/Databases/Vancouver, British Columbia.gpkg'
PARKS = gpd.read_feather('data/feather/cov_parks.feather')
REAL_TREES = gpd.read_feather('data/feather/nature_tree_canopy.feather')
land_use_gdf = gpd.read_feather('data/feather/mvan_landuse.feather')
diagonal_gdf = gpd.read_file('data/diagonal_tiles.geojson')


def gdf_from_memory(memory):
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
		getFillColor=[0, 0, 255, 100],
		get_line_color=[255, 255, 255, 255],
		auto_highlight=True,
		pickable=True,
	)
	if cell_id not in [l.id for l in deck.layers]:
		deck.layers.append(new_layer)
	else:
		deck.layers = [l for l in deck.layers if l.id != cell_id]
	return deck


# Define initial view state and create deck object
view_state = pdk.ViewState(latitude=49.254, longitude=-123.13, zoom=11, max_zoom=16, pitch=0, bearing=0)
r = pdk.Deck(
	layers=[],
	initial_view_state=view_state,
	map_style=pdk.map_styles.DARK_NO_LABELS,
)
mapbox_key = "pk.eyJ1IjoibmljaG9sYXNtYXJ0aW5vIiwiYSI6ImNrMjVhOGphOTAzZGUzbG8wNHJhdTZrMmYifQ.98uDMnGIvn1zrw4ZWUO35g"

tooltip = {
	"html": "F.S.R. {Max FSR} <br> Height {Max Height}"
}

html.Button(id='rotate', value='Rotate'),
right_bar_width = 380
layout = html.Div(
	children=[
		dcc.Store(id='memory'),
		Col(
			style={'height': '97vh', 'width': '100%', 'offset': 0, 'display': 'inline-block'},
			children=[

				# 3D Map
				Col(
					width={'size': 6, 'offset': 0},
					children=[
						html.Div(
							id="deck_div",
						),
					],
				),

				# Right panel
				Col(
					style={'width': f'{right_bar_width}px', 'height': f'90%', 'float': 'left',
					       'display': 'inline-block', 'font-family': 'roboto_light'},
					className='pretty_container',
					children=[
						dcc.Input(id='tiles_dir', value=GRID_FILE, style={"width": "100%"}),
						# dcc.Graph(id='sel_type_map'),
						html.Br(),
						html.Br(),
						# dl.Map(id=MAP_ID, center=[49.264, -123.115], zoom=16, children=[
						# 	dl.TileLayer(),  # Map tiles, defaults to OSM
						# 	dl.Polyline(id=POLYLINE_ID, positions=[dummy_pos]),
						# 	# Create a polyline, cannot be empty at the moment
						# 	dl.Polygon(id=POLYGON_ID, positions=[dummy_pos]),
						# 	# Create a polygon, cannot be empty at the moment
						# ], style={'width': '100%', 'height': '62%'}),
						html.Br(),
						html.Button('Select', id='select', style={"width": "37%"}),
						html.Button('Rotate', id='rotate', style={"width": "21%"}),
						html.Button('Filp H', id='flip_h', style={"width": "21%"}),
						html.Button('Flip V', id='flip_v', style={"width": "21%"}),
						html.Br(),
						dcc.Dropdown(id='type', options=[{'label': name, 'value': name} for name in types.keys()]),
					]
				),

				# Col(style={'height': f'680px', 'display': 'inline-block', 'float': 'right',}, children=[
				# 	# Layers menu
				# 	Row(
				# 		style={'width': f'{right_bar_width / 2}px', 'display': 'inline-block'},
				# 		children=[
				# 			html.Button(id='get_shp', value='SHP'),
				# 		]
				# 	),
				# ])

			]
		)
	]
)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.layout = layout


# STORES
def assign_callback(app_, out, event):
	@app_.callback(
		Output(component_id=f"memory", component_property="data"),
		[Input(component_id="deck", component_property=event)],
	)
	def dump_json(data):
		return data


assign_callback(app, "click-info", "clickInfo")


# LAYERS + DESIGN
@app.callback(
	Output(component_id="deck_div", component_property="children"),
	Output(component_id="area_by_lu", component_property="figure"),
	[
		Input(component_id='tiles_dir', component_property='value'),
		Input(component_id="type", component_property="value"),
		Input(component_id='rotate', component_property='n_clicks'),
		Input(component_id='flip_h', component_property='n_clicks'),
		Input(component_id='flip_v', component_property='n_clicks'),
		Input(component_id='select', component_property='n_clicks'),
		Input(component_id='memory', component_property='data'),
	])  # , State('input-on-submit', 'value'))
def main_callback(tiles_dir, sel_type, rotate, flip_h, flip_v, select, memory):
	stt = time.time()
	types_named = {t: i for i, t in enumerate(types)}

	# Test slash side
	if '/' in tiles_dir:
		splitter = '/'
	else:
		splitter = '\\'

	# Get filename
	file_name = tiles_dir.split(splitter)[-1:][0].split('.')[:-1][0]

	# # Test feather
	# if f'{file_name}.feather' not in os.listdir('data'):
	# 	grid_gdf = gpd.read_file(tiles_dir).to_crs(4326)
	# 	grid_gdf.to_feather(f'data/{file_name}.feather')
	# else:
	# 	grid_gdf = gpd.read_feather(f'data/{file_name}.feather')

	# Set index and random seed
	grid_gdf['id'] = grid_gdf.index
	np.random.seed(0)

	# Get callback context
	ctx = dash.callback_context
	if not ctx.triggered:
		button_id = 'No clicks yet'
	else:
		button_id = ctx.triggered[0]['prop_id'].split('.')[0]
	print(button_id)

	# Create cells deck layer
	if 'cells' not in [l.id for l in r.layers]:
		cells = pdk.Layer(
			id=f"cells",
			type="GeoJsonLayer",
			data=grid_gdf.to_crs(4326),
			stroked=True,
			wireframe=True,
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

	elif button_id == 'select':
		# Get selected cells
		for i in [l.id for l in r.layers]:
			if i.isnumeric():

				sel_cell = grid_gdf.loc[[int(i)], :].to_crs(26910)

				if sel_type == 'Typical Van SF':
					tp_gdf = gpd.read_file("/Users/nicholasmartino/Desktop/Sandbox Tiles/Typical_Van_SF.geojson")

				tp_uu = tp_gdf.to_crs(26910).unary_union
				tp_ctr = tp_uu.centroid
				cell_ctr = sel_cell.unary_union.centroid

				# Get difference between origin and destination geometry
				x_off = cell_ctr.x - tp_ctr.x
				y_off = cell_ctr.y - tp_ctr.y

				# Move one shape to a certain point
				tp_gdf = gpd.GeoDataFrame({
					'geometry': [translate(geom, x_off, y_off) for geom in tp_gdf['geometry']]}, crs=26910)

				r.layers.append(pdk.Layer(
					id=f"{sel_type}_{i}",
					type="GeoJsonLayer",
					data=tp_gdf.to_crs(4326),
					stroked=True,
					wireframe=True,
					getFillColor=[180, 180, 0, 255],
					get_line_color=[0, 0, 0],
					get_line_width=0.5,
					auto_highlight=True,
					pickable=False
				))

	# Run grids
	grid_gdf['id'] = grid_gdf.index
	grid_gdf['Type'] = grid_gdf['clus_gmm'].replace(TYPES)
	grid_gdf.loc[grid_gdf['Type'].isin(['Mid_High_Street', 'Moderate_Density', 'Dense_Nodal']), 'High St Type'] = 1
	STREETS['geometry'] = STREETS.buffer(5)
	grid_gdf.loc[gpd.overlay(grid_gdf, STREETS[STREETS['Category'] == 'Arterial'])['id'], 'Arterial'] = 1
	grid_gdf.loc[(grid_gdf['Arterial'] == 1) & (grid_gdf['High St Type'] == 1), 'High St'] = 1

	if memory is None:
		# Read/export all layers
		prefix = GRID_FILE.split('.')[0]
		all_layers_file = f"data/feather/{prefix}_all_layers.feather"
		if os.path.exists(all_layers_file):
			tiles = gpd.read_feather(all_layers_file)
		else:
			grid = Grid(grid_gdf, TILES, prefix=f"{prefix}_", land_use=land_use_gdf, diagonal_gdf=diagonal_gdf)
			grid.test_assign_subtypes()
			tiles = grid.test_place_tiles()
			tiles.to_feather(all_layers_file)

		prefix = f"{prefix.split('.')[0]}_"

		p_exists = os.path.exists(f"data/feather/{prefix}parcels.feather")
		b_exists = os.path.exists(f"data/feather/{prefix}buildings.feather")
		t_exists = os.path.exists(f"data/feather/{prefix}trees.feather")

		if p_exists and b_exists and t_exists:
			parcels = gpd.read_feather(f'data/feather/{prefix}parcels.feather')
			buildings = gpd.read_feather(f'data/feather/{prefix}buildings.feather')
			trees = gpd.read_feather(f'data/feather/{prefix}trees.feather')

		else:
			print(f"\n {prefix}")
			tiles = tiles.reset_index(drop=True)
			prcls = tiles[tiles['Type'] == 'prcls']
			bldgs = tiles[tiles['Type'] == 'bldgs']
			strts = tiles[tiles['Type'] == 'strts']
			blcks = tiles[tiles['Type'] == 'block']
			trees = tiles[tiles['Type'] == 'trees']

			scn = Scenario(parcels=prcls, buildings=bldgs, trees=trees, real_parks=PARKS, real_trees=REAL_TREES,
			               name=prefix)
			scn.parcels = scn.extract_parks()
			sb_trees = scn.extract_trees(directory='data/feather')
			sb_trees.to_feather(f"data/feather/{prefix}trees.feather")

			ind = Indicators(parcels=scn.parcels, buildings=bldgs, streets=strts, blocks=blcks)
			ind.test_indicators()
			ind.parcels.to_feather(f"data/feather/{prefix}parcels.feather")
			ind.buildings.loc[:, [c for c in ind.buildings.columns if c not in ['comm_units']]]. \
				to_feather(f"data/feather/{prefix}buildings.feather")
			ind.get_area_by_land_use().to_csv(f'data/feather/{prefix}land_use_area.csv')
			ind.get_floor_area_by_land_use().to_csv(f'data/feather/{prefix}land_use_floor_area.csv')
			ind.get_n_units_by_land_use().to_csv(f'data/feather/{prefix}land_use_n_units.csv')
			# stt = ind.get_street_length()
			# stt.loc[:, ["length", "geometry"]].to_file(f'{OUT_DIR}/{prefix}network_indicator.shp')
			# stt.loc[:, ["length"]].to_csv(f'{OUT_DIR}/{prefix}street_length.csv')
			# blocks = ind.get_block_area()
			# blocks.to_csv(f'{OUT_DIR}/{prefix}block_area.csv')
			# blocks.to_file(f'{OUT_DIR}/{prefix}blocks_indicator.shp')

			parcels = ind.parcels.copy()
			buildings = ind.buildings.copy()
			trees = sb_trees.copy()

			# all_tiles = pd.concat([ind.parcels, ind.buildings, sb_trees])
			# all_tiles.to_file(f'{OUT_DIR}/{prefix}all_tiles.geojson', driver='GeoJSON')
			# all_tiles.to_crs(4326).to_file(f'{OUT_DIR}/{prefix}all_tiles_4326.geojson', driver='GeoJSON')

			print(f"Total population: {ind.get_total_population()}")
			print(f"Total area: {ind.get_total_area()}")

		buildings['height'] = buildings['maxstories'] * 3

		# Create buildings and parcels layers
		parcels_pdk = pdk.Layer(
			id=f"parcels",
			type="GeoJsonLayer",
			data=parcels.loc[:, ['geometry']].to_crs(4326),
			getFillColor=[180, 180, 180, 255],
		)
		buildings_pdk = pdk.Layer(
			id=f"buildings",
			type="GeoJsonLayer",
			extruded=True,
			getElevation="height",
			data=buildings.loc[:, ['height', 'geometry']].to_crs(4326),
			getFillColor=[180, 180, 180, 255],
		)

		r.layers.append(parcels_pdk)
		r.layers.append(buildings_pdk)

		area_by_lu = None
	else:
		prefix = GRID_FILE.split('.')[0]
		prefix = f"{prefix.split('.')[0]}_"
		parcels = gpd.read_feather(f'data/feather/{prefix}parcels.feather')
		buildings = gpd.read_feather(f'data/feather/{prefix}buildings.feather')

		cells = gpd.GeoDataFrame(memory['object'], crs=4326).to_crs(26910)
		in_cell_bld = gpd.overlay(buildings, cells.loc[:, ['geometry']])
		in_cell_pcl = gpd.overlay(parcels, cells.loc[:, ['geometry']])

		ind = Indicators(parcels=in_cell_pcl, buildings=in_cell_bld)
		area_by_lu = px.bar(ind.get_area_by_land_use(), x='area', y='LANDUSE')
		ind.get_residential_units()
		ind.get_resident_count()

	# Create deck-gl object
	dgl = DeckGL(
		id="deck",
		data=r.to_json(),
		mapboxKey=mapbox_key,
		enableEvents=['click'],
		style={'width': '100%', 'float': 'left', 'display': 'inline-block'},
	)

	print(f"Callback: {round((time.time() - stt), 3)} seconds with {[l.id for l in r.layers]} layers")
	return dgl, area_by_lu


if __name__ == '__main__':
	app.run_server(debug=False, host='localhost', port=9000)
