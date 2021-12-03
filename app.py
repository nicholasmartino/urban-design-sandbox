import dash
from layout.layout import layout

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.layout = layout

server = app.server
app.config.suppress_callback_exceptions = True
