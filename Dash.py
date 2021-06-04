import time
import os
import dash
import dash_core_components as dcc
import dash_html_components as html
import geopandas as gpd
import numpy as np
import pydeck as pdk
from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate
from dash_bootstrap_components import Col, Row
from dash_deck import DeckGL
from shapely.geometry import Point, Polygon
import dash_leaflet as dl
from shapely.affinity import translate

alpha = 100
random_seed = 0
transparent = 'rgba(255,255,255,0)'

# Dash Leaflet parameters
MAP_ID = "map-id"
POLYLINE_ID = "polyline-id"
POLYGON_ID = "polygon-id"
dummy_pos = [0, 0]
dlatlon2 = 1e-6  # Controls tolerance of closing click

# Define zoning layer colors
types={
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


def gdf_from_memory(memory):
	memory['object']['geometry'] = [Polygon(coord) for coord in memory['object']['geometry']['coordinates']]
	gdf = gpd.GeoDataFrame(memory['object'], geometry='geometry', crs=4326)
	if 'id' not in gdf.columns:
		gdf['id'] = gdf.reset_index(drop=True).index
	return gdf

def pick_geometry(deck, gdf, cell_id):
	new_layer = pdk.Layer(
		id=cell_id,
		type="GeoJsonLayer",
		data=gdf,
		extruded=False,
		getFillColor=[0, 0, 255, 100],
		get_line_color=[255, 255, 255, 255],
		auto_highlight=True,
		pickable=True,
	)
	if cell_id not in [l.id for l in deck.layers]: deck.layers.append(new_layer)
	else: deck.layers = [l for l in deck.layers if l.id != cell_id]
	return deck


# Define initial view state and create deck object
view_state = pdk.ViewState(latitude=49.254, longitude=-123.13, zoom=11, max_zoom=16, pitch=0, bearing=0)
r = pdk.Deck(
	layers=[],
	initial_view_state=view_state,
	map_style=pdk.map_styles.DARK_NO_LABELS,
	mapbox_key="pk.eyJ1IjoibmljaG9sYXNtYXJ0aW5vIiwiYSI6ImNrMjVhOGphOTAzZGUzbG8wNHJhdTZrMmYifQ.98uDMnGIvn1zrw4ZWUO35g"
)


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
							style={'.bs-tooltip-left': {'top': '-12px', 'left': '-182px'}}
						),
					],
				),

				# Right panel
				Col(
					style={'width': f'{right_bar_width}px', 'height': f'90%', 'float': 'left',
						   'display': 'inline-block', 'font-family': 'roboto_light'},
					className='pretty_container',
					children=[
						dcc.Input(id='tiles_dir', value='/Volumes/Macintosh HD/Users/nicholasmartino/Google Drive/elementslab/CoV_gmm_r4_built.shp', style={"width": "100%"}),
						# dcc.Graph(id='sel_type_map'),
						html.Br(),
						html.Br(),
						dl.Map(id=MAP_ID, center=[49.264, -123.115], zoom=16, children=[
								dl.TileLayer(),  # Map tiles, defaults to OSM
								dl.Polyline(id=POLYLINE_ID, positions=[dummy_pos]),  # Create a polyline, cannot be empty at the moment
								dl.Polygon(id=POLYGON_ID, positions=[dummy_pos]),  # Create a polygon, cannot be empty at the moment
							], style={'width': '100%', 'height': '62%'}),
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


### STORES ###
def assign_callback(app_, out, event):
	@app_.callback(
		Output(component_id=f"memory", component_property="data"),
		[Input(component_id="deck", component_property=event)],
	)
	def dump_json(data):
		return data
assign_callback(app, "click-info", "clickInfo")


@app.callback([Output(POLYLINE_ID, "positions"), Output(POLYGON_ID, "positions")],
			  [Input(MAP_ID, "click_lat_lng")],
			  [State(POLYLINE_ID, "positions")])
def update_polyline_and_polygon(click_lat_lng, positions):
	if click_lat_lng is None or positions is None:
		raise PreventUpdate()
	# On first click, reset the polyline.
	if len(positions) == 1 and positions[0] == dummy_pos:
		return [click_lat_lng], [dummy_pos]
	# If the click is close to the first point, close the polygon.
	dist2 = (positions[0][0] - click_lat_lng[0]) ** 2 + (positions[0][1] - click_lat_lng[1]) ** 2
	if dist2 < dlatlon2:
		return [dummy_pos], positions
	# Otherwise, append the click position.
	positions.append(click_lat_lng)
	return positions, [dummy_pos]

	Polygon([Point(c) for c in positions])


### LAYERS + DESIGN ###
@app.callback(
	Output(component_id="deck_div", component_property="children"),
	[
	 Input(component_id='tiles_dir', component_property='value'),
	 Input(component_id="type", component_property="value"),
	 Input(component_id='rotate', component_property='n_clicks'),
	 Input(component_id='flip_h', component_property='n_clicks'),
	 Input(component_id='flip_v', component_property='n_clicks'),
	 Input(component_id='select', component_property='n_clicks'),
	 Input(component_id='memory', component_property='data'),
	])#, State('input-on-submit', 'value'))
def main_callback(tiles_dir, sel_type, rotate, flip_h, flip_v, select, memory):
	stt = time.time()
	types_named = {t: i for i, t in enumerate(types)}

	# Test slash side
	if '/' in tiles_dir: splitter = '/'
	else: splitter = '\\'

	# Get filename
	file_name = tiles_dir.split(splitter)[-1:][0].split('.')[:-1][0]

	# Test feather
	if f'{file_name}.feather' not in os.listdir('data'):
		grid_gdf = gpd.read_file(tiles_dir).to_crs(4326)
		grid_gdf.to_feather(f'data/{file_name}.feather')
	else:
		grid_gdf = gpd.read_feather(f'data/{file_name}.feather')

	# Set index and random seed
	grid_gdf['id'] = grid_gdf.index
	np.random.seed(0)

	# Get callback context
	ctx = dash.callback_context
	if not ctx.triggered: button_id = 'No clicks yet'
	else: button_id = ctx.triggered[0]['prop_id'].split('.')[0]
	print(button_id)

	# Create cells deck layer
	if 'cells' not in [l.id for l in r.layers]:
		cells = pdk.Layer(
			id=f"cells",
			type="GeoJsonLayer",
			data=grid_gdf,
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

				if sel_type == 'Typical Van SF': tp_gdf = gpd.read_file("/Users/nicholasmartino/Desktop/Sandbox Tiles/Typical_Van_SF.geojson")

				tp_uu = tp_gdf.to_crs(26910).unary_union
				tp_ctr = tp_uu.centroid
				cell_ctr = sel_cell.unary_union.centroid

				# Get difference between origin and destination geometry
				x_off = cell_ctr.x - tp_ctr.x
				y_off = cell_ctr.y - tp_ctr.y

				# Move one shape to a certain point
				tp_gdf = gpd.GeoDataFrame({'geometry': [translate(geom, x_off, y_off) for geom in tp_gdf['geometry']]}, crs=26910)

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

	# Create deck-gl object
	dgl = DeckGL(
		id="deck",
		data=r.to_json(),
		mapboxKey=r.mapbox_key,
		enableEvents=['click'],
		style={'width': '100%', 'float': 'left', 'display': 'inline-block'},
	)
	print(f"Callback: {round((time.time() - stt), 3)} seconds with {[l.id for l in r.layers]} layers")
	return dgl



if __name__ == '__main__':
	app.run_server(debug=False, host='localhost', port=9000)
