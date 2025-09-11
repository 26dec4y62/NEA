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
        QtWidgets.QMessageBox.information(self, 'Works', f'{name} recognised.')
        self.appdata.add_history(name)
        self.appdata.fetch_molecule(name)

def main():
    app = QtWidgets.QApplication(sys.argv)
    win = MoleculeViewer()
    win.show()
    sys.exit(app.exec())

# if __name__ == '__main__':
main()