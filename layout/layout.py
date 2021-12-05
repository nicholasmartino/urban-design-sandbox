<<<<<<< HEAD
import dash_bootstrap_components as dbc
=======
from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
>>>>>>> ef39bd4ecab44754cfa3e028344731276ce4a7af
import dash_extensions as de
import dash_core_components as dcc
import dash_html_components as html

RIGHT_BAR_WIDTH = 380
# GRID_FILE = "fishned_CoV_gmm_joined.geojson"
ALPHA = 100
types = {
	'Open_Low_Density': 2,
	'Mid_High_Street': 3,
	'Coarse_Grain': 4,
	'Moderate_Density': 5,
	'Treed_Large_Homes': 8,
	'Dense_Nodal': 9,
	'Typical_Van_SF': 10,
	'Typical_Van_SF_West': 11,
	'Open Green Space': 12,
}
TEXT_INDICATOR = {'opacity': '1.0', 'margin-left': '30px'}

layout = html.Div(
	children=[
		dcc.Store(id='memory'),
		dbc.Col(
			style={'height': '99vh', 'width': '100%', 'offset': 0, 'display': 'inline-block'},
			children=[

				# 3D Map
				dbc.Col(
					width={'size': 6, 'offset': 0},
					children=[
						html.Div(id="deck_div"),
						dcc.Loading(
							id="loading-1",
							type="default",
							children=html.Div(id="loading-output-1"),
							style={"position": "fixed", 'top': '50%', 'left': '50%'}
						)
					]
				),

				# Left panel
				dbc.Col(
					style={
						'width': f'{RIGHT_BAR_WIDTH}px', 'height': f'95%', 'float': 'left',
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
						html.Div(id="selected_type", style={"width": "100%", 'text-align': 'center'}),
						dcc.Graph(id='cell_type_icon', style={'height': '20%'}),
						dcc.Graph(id='area_by_lu', style={'height': '40%'}),
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

				# Right download button
				dbc.Col(
					style={
						'width': f'190px', 'height': f'40px', 'float': 'right',
						'display': 'inline-block', 'font-family': 'Avenir'
					},
					className='pretty_container',
					children=[
						html.Button("Download Layers", id="btn-download", style={
							"float": "right", 'display': 'inline-block'
						}),
						de.Download(id="download")
					]
				)
			]
		)
	]
)
