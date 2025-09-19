import json
import requests
from urllib.parse import quote

FILE = 'molecules.json'

class appdata:
    def __init__(self):
        self.names = []
        self.properties = []
        self.load_data()

    def load_data(self):
        with open(FILE, 'r') as f:
            data = json.load(f)
        self.names = data.get('History', [])
        self.properties = data.get('Properties', [])

    def save_data(self):
        with open(FILE, 'w') as f:
            json.dump({'History': self.names, 'Properties': self.properties}, f, indent=2)
        return

    def add_history(self, name):
        if name not in self.names:
            self.names.append(name)
            self.save_data()

    def fetch_molecule(self, name):
        name_enc = quote(name)
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name_enc}property/ConnectivitySMILES,IUPACName/JSON"
        )
        r = requests.get(url, timeout=15)
        r.raise_for_status()
        js = r.json()
        props_list = js.get("PropertyTable", {}).get("Properties", [])
        print(props_list)
        if not props_list:
            raise ValueError(f'No properties for {name}')
        prop = props_list[0]
        self.properties = [prop]
        self.save_data()
        return prop