import sys
from PyQt5 import QtWidgets
import vtk
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from backend import appdata
from vtk_creation import MoleculeRenderer

class Window(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecule Viewer")
        self.resize(1000, 800)
        self.appdata = appdata()

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)

        # Input section
        input_layout = QtWidgets.QHBoxLayout()
        self.entry = QtWidgets.QLineEdit()
        self.entry.setPlaceholderText("Enter molecule name e.g., methanol")
        self.entry.returnPressed.connect(self.search)
        btn = QtWidgets.QPushButton("Search & Render")
        btn.clicked.connect(self.search)
        input_layout.addWidget(self.entry)
        input_layout.addWidget(btn)
        layout.addLayout(input_layout)

        # Info section
        info_layout = QtWidgets.QVBoxLayout()
        self.iupac_name = QtWidgets.QLabel("IUPAC Name: —")
        self.structural_formula = QtWidgets.QLabel("Molecular Formula: —")
        self.atom_count = QtWidgets.QLabel("Atoms: —")
        info_layout.addWidget(self.iupac_name)
        info_layout.addWidget(self.structural_formula)
        info_layout.addWidget(self.atom_count)
        layout.addLayout(info_layout)

        # VTK widget for 3D rendering
        self.vtk_widget = QVTKRenderWindowInteractor(self)
        layout.addWidget(self.vtk_widget)
        
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(0.1, 0.1, 0.1)
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)
        self.iren = self.vtk_widget.GetRenderWindow().GetInteractor()
        
        # Set up interactor style for better control
        style = vtk.vtkInteractorStyleTrackballCamera()
        self.iren.SetInteractorStyle(style)
        self.iren.Initialize()

        # Initialize molecule renderer
        self.molecule_renderer = MoleculeRenderer(self.renderer)

    def search(self):
        molecule_name = self.entry.text().strip()
        if not molecule_name:
            QtWidgets.QMessageBox.warning(self, "Input Error", "Please enter a molecule name")
            return
        
        try:
            # Get molecule data from backend
            molecule_data = self.appdata.get_molecule(molecule_name)
            
            # Update labels
            self.iupac_name.setText(f"IUPAC Name: {molecule_data.get("iupac_name", "—")}")
            self.structural_formula.setText(f"Structural Formula: {molecule_data.get("structural_formula", "—")}")
            
            # Render molecule using vtk_creation module
            success = self.molecule_renderer.render_molecule(molecule_data)
            
            if success:
                atom_count = self.molecule_renderer.get_atom_count()
                bond_count = self.molecule_renderer.get_bond_count()
                self.atom_count.setText(f"Atoms: {atom_count} | Bonds: {bond_count}")
                self.update_display()
                self.statusBar().showMessage(f"Successfully loaded {molecule_name}")
            else:
                QtWidgets.QMessageBox.warning(self, "No Data", 
                    "No atom data available for this molecule")
                self.statusBar().showMessage("No data available")
            
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", 
                f"Failed to load molecule:\n{str(e)}")
    
    def update_display(self):
        self.vtk_widget.GetRenderWindow().Render()

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec())