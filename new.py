import sys
import os
from PyQt5 import QtWidgets
import vtk
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import json
import requests
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

FILE = 'molecules.json'

class backend:
    def load_data(self):
        if os.path.exists(FILE):
            with open(FILE, 'r') as f:
                self.data = json.load(f)
    
    def save_data(self):
        with open(FILE, 'w') as f:
            json.dump(self.data, f, indent=2)

class MoleculeViewer(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecule Viewer")
        self.resize(800, 600)
        self.backend = backend()

        # Central widget
        self.frame = QtWidgets.QFrame()
        self.layout = QtWidgets.QVBoxLayout()
        self.frame.setLayout(self.layout)
        self.setCentralWidget(self.frame)

        # VTK Widget
        self.vtk_widget = QVTKRenderWindowInteractor(self)
        self.layout.addWidget(self.vtk_widget)

        # Renderer
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(1, 1, 1)
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        self.iren.Initialize()

        # Interactor
        self.interactor = self.vtk_widget.GetRenderWindow().GetInteractor()
        self.interactor.Initialize()

        # Load molecule data
        self.backend.load_data()

def main():
    app = QtWidgets.QApplication(sys.argv)
    win = MoleculeViewer()
    win.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()