# Imports
from pyexpat import expat_CAPI
import sys
import math
import requests
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import vtk
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtWidgets
from urllib.parse import quote

# Atom colors & radii
ATOM_COLOURS = {
    "H": (1.0, 1.0, 1.0),
    "C": (0.2, 0.2, 0.2),
    "N": (0.0, 0.0, 1.0),
    "O": (1.0, 0.0, 0.0),
    "F": (0.0, 0.8, 0.0),
    "Cl": (0.0, 0.8, 0.0),
    "Br": (0.6, 0.2, 0.2),
    "I": (0.4, 0.0, 0.6),
    "S": (1.0, 0.8, 0.2),
    "P": (1.0, 0.5, 0.0),
}
ATOM_RADII = {
    "H": 0.25,
    "C": 0.35,
    "N": 0.33,
    "O": 0.33,
    "F": 0.32,
    "Cl": 0.38,
    "Br": 0.40,
    "I": 0.45,
    "S": 0.38,
    "P": 0.40
}

def color_for_symbol(sym):
    return ATOM_COLOURS.get(sym, (0.7, 0.7, 0.7))

def radius_for_symbol(sym):
    return ATOM_RADII.get(sym, 0.35)

def fetch_pubchem_by_name(name: str):
    # Fetch SMILES and IUPACName from PubChem by molecule name
    name_enc = quote(name)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name_enc}/property/CanonicalSMILES,IUPACName/JSON"
    
    r = requests.get(url, timeout=10)
    
    r.raise_for_status()
    js = r.json()
    
    props_list = js.get("PropertyTable", {}).get("Properties", [])
    if not props_list:
        raise ValueError(f"No properties returned by PubChem for '{name}'")
    
    props = props_list[0]
    smiles = props.get("CanonicalSMILES") or props.get("IsomericSMILES") or props.get("ConnectivitySMILES")
    iupac = props.get("IUPACName") or "(IUPAC unavailable)"
    
    if not smiles:
        raise ValueError(f"No SMILES returned by PubChem for '{name}'")
    
    return {"smiles": smiles, "iupac": iupac}

# RDKit 3D molecule
def mol_from_smiles_3d(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES.")
    mol = Chem.AddHs(mol)
    status = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if status != 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    return mol


# VTK rendering
def add_molecule_to_renderer(mol, renderer):
    renderer.RemoveAllViewProps()
    conf = mol.GetConformer()

    # Parsing atoms
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        sym = atom.GetSymbol()

        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(p.x, p.y, p.z)
        sphere.SetRadius(radius_for_symbol(sym))
        sphere.SetThetaResolution(24)
        sphere.SetPhiResolution(24)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(*color_for_symbol(sym))
        renderer.AddActor(actor)

    # Making bonds
    for bond in mol.GetBonds():
        start = conf.GetAtomPosition(bond.GetBeginAtomIdx())
        end = conf.GetAtomPosition(bond.GetEndAtomIdx())
        v = np.array([end.x - start.x, end.y - start.y, end.z - start.z])
        length = np.linalg.norm(v)
        if length == 0: continue

        # Create cylinder
        cyl = vtk.vtkCylinderSource()
        cyl.SetRadius(0.12)
        cyl.SetResolution(20)
        cyl.SetHeight(1)

        # Transform cylinder
        transform = vtk.vtkTransform()
        mid = [(start.x + end.x)/2, (start.y + end.y)/2, (start.z + end.z)/2]
        transform.Translate(*mid)

        # Alignment
        v_norm = v / length
        y_axis = np.array([0, 1, 0])
        axis = np.cross(y_axis, v_norm)
        angle = np.degrees(np.arccos(np.clip(np.dot(y_axis, v_norm), -1.0, 1.0)))
        if np.linalg.norm(axis) > 1e-6:
            transform.RotateWXYZ(angle, *axis)

        transform.Scale(1, length, 1)

        tf = vtk.vtkTransformPolyDataFilter()
        tf.SetTransform(transform)
        tf.SetInputConnection(cyl.GetOutputPort())

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(tf.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0.3, 0.3, 0.3)
        renderer.AddActor(actor)

    renderer.ResetCamera()
    renderer.GetActiveCamera().Zoom(1.3)


# PyQt GUI
class MoleculeViewer(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("3D Molecule Display")
        self.resize(1000, 700)

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)

        # Input row
        input_layout = QtWidgets.QHBoxLayout()
        self.entry = QtWidgets.QLineEdit()
        self.entry.setPlaceholderText("Enter molecule name")
        btn = QtWidgets.QPushButton("Search & Render")
        btn.clicked.connect(self.on_search)
        input_layout.addWidget(self.entry)
        input_layout.addWidget(btn)
        layout.addLayout(input_layout)

        # Labels
        self.iupac_label = QtWidgets.QLabel("IUPAC: —")
        self.smiles_label = QtWidgets.QLabel("SMILES: —")
        layout.addWidget(self.iupac_label)
        layout.addWidget(self.smiles_label)

        # VTK widget
        self.vtkWidget = QVTKRenderWindowInteractor(self)
        layout.addWidget(self.vtkWidget)

        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(1, 1, 1)
        self.vtkWidget.GetRenderWindow().AddRenderer(self.renderer)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        self.iren.Initialize()

    def on_search(self):
        name = self.entry.text().strip()
        if not name:
            QtWidgets.QMessageBox.warning(self, "Input error", "Please enter a molecule name.")
            return
        try:
            data = fetch_pubchem_by_name(name)
            smiles, iupac = data["smiles"], data["iupac"]
            self.iupac_label.setText(f"IUPAC: {iupac}")
            self.smiles_label.setText(f"SMILES: {smiles}")

            mol = mol_from_smiles_3d(smiles)
            add_molecule_to_renderer(mol, self.renderer)
            self.vtkWidget.GetRenderWindow().Render()

        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", str(e))

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    win = MoleculeViewer()
    win.show()
    sys.exit(app.exec())