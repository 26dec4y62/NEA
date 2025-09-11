import os
import json
import requests
from urllib.parse import quote

FILE = 'molecules.json'

class appdata:
    def __init__(self):
        self.names = []
        self.properties = []

    def load_data(self):
        if os.path.exists(FILE):
            with open(FILE, 'r') as f:
                json.load(f)
    
    def save_data(self):
        with open(FILE, 'w') as f:
            json.dump({'History': self.names, 'Properties': self.properties}, f, indent=2)

    def add_history(self, name):
        if name and name not in self.names:
            self.names.append(name)
            self.save_data()

    def fetch_molecule(self, name):
        name_enc = quote(name)
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name_enc}/property/MolecularFormula/JSON"
        )
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        js = r.json()
        self.properties = js.get("PropertyTable", {}).get("Properties", [])
        if not self.properties:
            raise ValueError(f'No properties for {name}')
        
        props = self.properties[0]
        molecular = (props.get('MolecularFormula', ''))
        iupac_name = props.get('IUPACName', '')
        self.properties.append([molecular, iupac_name])
        self.save_data()
