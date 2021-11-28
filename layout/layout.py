import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html

RIGHT_BAR_WIDTH = 380
GRID_FILE = "fishned_CoV_gmm_joined.geojson"
ALPHA = 100
types = {
	'Typical Van SF': [255, 100, 100, ALPHA],
	'2': [255, 100, 255, ALPHA],
	'3': [100, 100, 255, ALPHA],
	'4': [255, 255, 100, ALPHA],
	'5': [255, 127, 100, ALPHA],
	'6': [255, 127, 100, ALPHA],
	'7': [255, 127, 100, ALPHA],
	'8': [255, 127, 100, ALPHA],
	'9': [255, 127, 100, ALPHA],
	'10': [255, 127, 100, ALPHA],
	'11': [255, 127, 100, ALPHA],
}
TEXT_INDICATOR = {'opacity': '1.0', 'margin-left': '30px'}

layout = html.Div(
	children=[
		dcc.Store(id='memory'),
		dbc.Col(
			style={'height': '100vh', 'width': '100%', 'offset': 0, 'display': 'inline-block'},
			children=[

				# 3D Map
				dbc.Col(
					width={'size': 6, 'offset': 0},
					children=[html.Div(id="deck_div")]
				),

				# Right panel
				dbc.Col(
					style={
						'width': f'{RIGHT_BAR_WIDTH}px', 'height': f'90%', 'float': 'left',
						'display': 'inline-block', 'font-family': 'Avenir'
					},
					className='pretty_container',
					children=[
						dcc.Dropdown(id='type', options=[{'label': name, 'value': name} for name in types.keys()]),
						dbc.Row([
							html.Button('Select', id='select', style={"width": "37%"}),
							html.Button('Rotate', id='rotate', style={"width": "21%"}),
							html.Button('Filp H', id='flip_h', style={"width": "21%"}),
							html.Button('Flip V', id='flip_v', style={"width": "21%"}),
						]),
						html.Br(),
						dcc.Graph(id='cell_type_icon', style={'height': '31%'}),
						dcc.Graph(id='area_by_lu', style={'height': '31%'}),
						dbc.Row([
							html.Div(id='total_units', style=TEXT_INDICATOR)
						]),
						dbc.Row([
							html.Div(id='total_population', style=TEXT_INDICATOR)
						]),
						html.Br(),
						html.H5("Upload grid file"),
						dcc.Upload(
							id="upload",
							children=html.Div([
								'Drag and Drop or ',
								html.A('Select Files')
							]),
							style={
								'width': '95%',
								'height': '50%',
								'lineHeight': '60px',
								'borderWidth': '1px',
								'borderStyle': 'dashed',
								'borderRadius': '5px',
								'textAlign': 'center',
								'margin': '10px'
							},
						),
						html.Small("Accepted formats: GeoJSON, Shapefile. Must contain 'clus_gmm' and 'Subtype' columns")
					]
				),
			]
		)
	]
)
