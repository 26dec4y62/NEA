import sys
from PyQt5 import QtWidgets, QtGui
import vtk
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
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
        try:
            self.key_window.hide()
        except Exception:
            pass
        name = self.entry.text().strip()
        if not name:
            QtWidgets.QMessageBox.warning(self, 'Input error', 'Please enter a molecule name.')
            return
        QtWidgets.QMessageBox.information(self, 'Check', f'{name} recognised.')
        self.display_key()
        self.appdata.add_history(name)
        try:
            data = self.appdata.fetch_molecule(name)
            print(data)
            self.iupac, self.structural = data['IUPACName'], Molecule.structural_formula(data['SMILES'])
        except Exception:
            QtWidgets.QMessageBox.critical(self, "Error", f'Could not fetch data for {name}.')
            return
        # Molecule.parse_molecule_to_dict(data['SMILES'])
        self.update_display(self.iupac, self.structural)

    def update_display(self, iupac, structural):
        self.iupac_name.setText(f'IUPAC Name: {iupac}')
        self.structural_formula.setText(f'Structural Formula: {structural}')

    def display_key(self):
        # Create a list widget to show atom colors
        self.key_window = QtWidgets.QDialog(self)
        self.key_window.setWindowTitle('Atom Color Key')
        self.key_window.resize(300, 400)
        
        layout = QtWidgets.QVBoxLayout(self.key_window)
        key_list = QtWidgets.QListWidget()
        
        # Add each atom with its color information
        lst = [char for char in set(self.structural) if char.isalpha()]
        for symbol, color in Molecule.ATOM_COLOURS.items():
            rgb_text = f"RGB({color[0]:.1f}, {color[1]:.1f}, {color[2]:.1f})"
            item_text = f"{symbol} - {rgb_text}"
            item = QtWidgets.QListWidgetItem(item_text)
            r, g, b = [int(c * 255) if symbol != "H" else int(0.7*255) for c in color] # Light grey for H
            item.setBackground(QtGui.QColor(r, g, b))
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