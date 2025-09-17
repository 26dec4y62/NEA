import sys
from PyQt5 import QtWidgets
import vtk
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
# import numpy as np
# from rdkit import Chem
# from rdkit.Chem import AllChem
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
        self.show_molecule(name)
        self.display_key()
        self.appdata.add_history(name)
        try:
            data = self.appdata.fetch_molecule(name)
            iupac, structural = data['IUPACName'], Molecule.structural_formula(data['MolecularFormula'])
            print(iupac, structural)
            self.update_display(iupac, structural)
        except Exception:
            QtWidgets.QMessageBox.critical(self, "Error", f'Could not fetch data for {name}.')
            return

    def show_molecule(self, smiles):
        return

    def update_display(self, iupac, structural):
        self.iupac_name.setText(f'IUPAC Name: {iupac}')
        self.structural_formula.setText(f'Structural Formula: {structural}')
        return

    def display_key(self):
        key = QtWidgets.QMdiSubWindow()
        key.setWindowTitle('Key')
        key.resize(300, 200)
        # self.key_list = QtWidgets.QListWidget()
        # key.setWidget(self.key_list)   
        # for symbol in Molecule.ATOM_COLOURS.keys():
        #     self.key_list.addItem(symbol)
        key.show()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    win = MoleculeViewer()
    win.show()
    sys.exit(app.exec())