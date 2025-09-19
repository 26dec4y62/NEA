import sys
from PyQt5 import QtWidgets, QtGui, QtCore
import vtk
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import numpy as np
import math
# Connecting files
from backend import appdata
from chem import Molecule

class MoleculeViewer(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Molecule Viewer')
        self.resize(800, 600)
        self.appdata = appdata()

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)

        # Input row
        input_layout = QtWidgets.QHBoxLayout()
        self.entry = QtWidgets.QLineEdit()
        self.entry.setPlaceholderText('Enter molecule name (e.g. methanol)')
        btn = QtWidgets.QPushButton('Search & Render')
        btn.clicked.connect(self.on_search)
        input_layout.addWidget(self.entry)
        input_layout.addWidget(btn)
        layout.addLayout(input_layout)

        # Labels
        self.iupac_name = QtWidgets.QLabel('IUPAC Name: —')
        self.structural_formula = QtWidgets.QLabel('Structural Formula: —')
        layout.addWidget(self.iupac_name)
        layout.addWidget(self.structural_formula)

        # VTK Widget
        self.vtk_widget = QVTKRenderWindowInteractor(self)
        layout.addWidget(self.vtk_widget)

        # Renderer
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(1, 1, 1)
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)
        self.iren = self.vtk_widget.GetRenderWindow().GetInteractor()
        self.iren.Initialize()  

    def on_search(self):
        name = self.entry.text().strip()
        if not name:
            QtWidgets.QMessageBox.warning(self, 'Input error', 'Please enter a molecule name.')
            return
        QtWidgets.QMessageBox.information(self, 'Check', f'{name} recognised.')
        try:
            data = self.appdata.fetch_molecule(name) # fetch is returning error
            iupac, structural = data['IUPACName'], Molecule.structural_formula(data['ConnectivitySMILES'])
            self.update_display(iupac, structural)
            self.display_key()
            self.appdata.add_history(name)
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", f'Could not fetch data for {name}.')
            return

    def _extract_atoms_and_bonds_rdkit(self):
        if self.mol is None:
            self.atoms = []
            self.bonds = []
            return

        conf = self.mol.GetConformer()
        self.atoms = []
        for atom in self.mol.GetAtoms():
            idx = atom.GetIdx()
            pos = conf.GetAtomPosition(idx)
            self.atoms.append({
                'symbol': atom.GetSymbol(),
                'position': (pos.x, pos.y, pos.z),
                'index': idx
            })

        self.bonds = []
        for bond in self.mol.GetBonds():
            self.bonds.append({
                'a1': bond.GetBeginAtomIdx(),
                'a2': bond.GetEndAtomIdx(),
                'order': int(bond.GetBondTypeAsDouble())
            })

    # VTK rendering
    def create_atom_actors(self):
        actors = []
        for a in self.atoms:
            x, y, z = a['position']
            sym = a['symbol']
            r = self.ATOM_RADII.get(sym, 0.35)

            sphere = vtk.vtkSphereSource()
            sphere.SetCenter(x, y, z)
            sphere.SetRadius(r)
            sphere.SetThetaResolution(24)
            sphere.SetPhiResolution(24)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(sphere.GetOutputPort())

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            color = self.ATOM_COLOURS.get(sym, (0.7, 0.7, 0.7))
            actor.GetProperty().SetColor(*color)
            actors.append(actor)
        return actors

    def create_bond_actors(self):
        actors = []
        for b in self.bonds:
            a1 = self.atoms[b['a1']]
            a2 = self.atoms[b['a2']]
            start = np.array(a1['position'])
            end = np.array(a2['position'])
            v = end - start
            length = np.linalg.norm(v)
            if length == 0:
                continue

            cyl = vtk.vtkCylinderSource()
            cyl.SetRadius(0.12)
            cyl.SetResolution(24)
            cyl.SetHeight(1.0)

            # transform: translate to midpoint, rotate and scale
            transform = vtk.vtkTransform()
            mid = (start + end) / 2.0
            transform.Translate(*mid)

            # compute rotation
            y_axis = np.array([0.0, 1.0, 0.0])
            v_dir = v / length
            axis = np.cross(y_axis, v_dir)
            axis_norm = np.linalg.norm(axis)
            if axis_norm < 1e-8:
                angle = 0.0
                axis = (1.0, 0.0, 0.0)
            else:
                angle = math.degrees(math.acos(np.clip(np.dot(y_axis, v_dir), -1.0, 1.0)))
                axis = axis / axis_norm

            if angle != 0.0:
                transform.RotateWXYZ(angle, axis[0], axis[1], axis[2])

            transform.Scale(1.0, length, 1.0)

            tf = vtk.vtkTransformPolyDataFilter()
            tf.SetTransform(transform)
            tf.SetInputConnection(cyl.GetOutputPort())
            tf.Update()

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(tf.GetOutputPort())
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(0.35, 0.35, 0.35)
            actors.append(actor)
        return actors

    def update_display(self, iupac, structural):
        self.iupac_name.setText(f'IUPAC Name: {iupac}')
        self.structural_formula.setText(f'Structural Formula: {structural}')

    def display_key(self):
        self.key_window = QtWidgets.QDialog(self)
        self.key_window.setWindowTitle('Atom Color Key')
        self.key_window.resize(300, 400)
        
        layout = QtWidgets.QVBoxLayout(self.key_window)
        
        # Create a list widget to show atom colors
        key_list = QtWidgets.QListWidget()
        
        # Add each atom with its color information
        for symbol, color in Molecule.ATOM_COLOURS.items():
            rgb_text = f"RGB({color[0]:.1f}, {color[1]:.1f}, {color[2]:.1f})"
            item_text = f"{symbol} - {rgb_text}"
            item = QtWidgets.QListWidgetItem(item_text)
            r, g, b = [int(c * 255) for c in color]
            item.setBackground(QtWidgets.QColor(r, g, b))
            
            # Set text color to contrast with background
            if sum(color) < 1.5:  # Dark background
                item.setForeground(QtWidgets.QColor(255, 255, 255))  # White text
            else:  # Light background
                item.setForeground(QtWidgets.QColor(0, 0, 0))  # Black text
            
            key_list.addItem(item)
        layout.addWidget(key_list)
        close_btn = QtWidgets.QPushButton('Close')
        close_btn.clicked.connect(self.key_window.close)
        layout.addWidget(close_btn)
        self.key_window.show()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    win = MoleculeViewer()
    win.show()
    sys.exit(app.exec())