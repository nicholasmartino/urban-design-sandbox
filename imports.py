import requests


def import_github(user='nicholasmartino', repository='morphology', branch='master', file='main.py'):
    r = requests.get(f"https://raw.githubusercontent.com/{user}/{repository}/{branch}/{file}")
    with open(f"models/{file}", "wb") as f:
        f.write(r.content)
    return


import_github(repository='morphology', file='ShapeTools.py')
import_github(repository='city', file='Network.py')
