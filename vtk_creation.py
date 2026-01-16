import vtk
import math

class Actor():
    def __init__(self, actor=None):
        self.actor = actor if actor is not None else vtk.vtkActor()
    
    def set_color(self, r, g, b):
        self.actor.GetProperty().SetColor(r, g, b)
        return 
    
    def set_opacity(self, opacity):
        self.actor.GetProperty().SetOpacity(opacity)
        return
    
    def set_position(self, x, y, z):
        self.actor.SetPosition(x, y, z)
        return 
    
    def set_scale(self, sx, sy, sz):
        self.actor.SetScale(sx, sy, sz)
        return 
    
    def get_actor(self):
        return self.actor


class Atom(Actor):
    def __init__(self, element, position, radius=0.3):
        super().__init__()
        self.element = element
        self.position = position
        self.radius = radius
        
        # Create sphere geometry
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(position[0], position[1], position[2])
        sphere.SetRadius(radius)
        sphere.SetPhiResolution(20)
        sphere.SetThetaResolution(20)
        
        # Create mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
        
        # Set up actor
        self.actor.SetMapper(mapper)
        
        # Set colour based on element (CPK coloring)
        self._set_element_colour()
    
    def _set_element_colour(self):
        colours = {
            'H': (1.0, 1.0, 1.0),    # White
            'C': (0.5, 0.5, 0.5),    # Gray
            'N': (0.0, 0.0, 1.0),    # Blue
            'O': (1.0, 0.0, 0.0),    # Red
            'F': (0.0, 1.0, 0.0),    # Green
            'Cl': (0.0, 1.0, 0.0),   # Green
            'Br': (0.6, 0.2, 0.0),   # Brown
            'I': (0.5, 0.0, 0.5),    # Purple
            'S': (1.0, 1.0, 0.0),    # Yellow
            'P': (1.0, 0.5, 0.0),    # Orange
        }
        colour = colours.get(self.element, (1.0, 0.0, 1.0))  # Magenta for unknown
        self.set_color(*colour)


class Bond(Actor):
    def __init__(self, atom1_pos, atom2_pos, radius=0.1):
        super().__init__()
        
        # Create cylinder
        cylinder = vtk.vtkCylinderSource()
        cylinder.SetRadius(radius)
        cylinder.SetHeight(1.0)
        cylinder.SetResolution(20)
        
        # Calculate bond vector and length using plain Python
        v1_x, v1_y, v1_z = atom1_pos
        v2_x, v2_y, v2_z = atom2_pos
        
        bond_x = v2_x - v1_x
        bond_y = v2_y - v1_y
        bond_z = v2_z - v1_z
        
        bond_length = math.sqrt(bond_x**2 + bond_y**2 + bond_z**2)
        center_x = (v1_x + v2_x) / 2
        center_y = (v1_y + v2_y) / 2
        center_z = (v1_z + v2_z) / 2
        
        # Create transform to position and orient the cylinder
        transform = vtk.vtkTransform()
        transform.Translate(center_x, center_y, center_z)
        
        # Rotate to align with bond
        bond_unit_x = bond_x / bond_length
        bond_unit_y = bond_y / bond_length
        bond_unit_z = bond_z / bond_length
        
        # Cylinder's default axis is (0, 1, 0)
        cylinder_axis = (0, 1, 0)
        
        # Cross product for rotation axis
        rot_x = cylinder_axis[1] * bond_unit_z - cylinder_axis[2] * bond_unit_y
        rot_y = cylinder_axis[2] * bond_unit_x - cylinder_axis[0] * bond_unit_z
        rot_z = cylinder_axis[0] * bond_unit_y - cylinder_axis[1] * bond_unit_x
        
        rot_length = math.sqrt(rot_x**2 + rot_y**2 + rot_z**2)
        
        # Dot product for rotation angle
        dot_product = cylinder_axis[0] * bond_unit_x + cylinder_axis[1] * bond_unit_y + cylinder_axis[2] * bond_unit_z
        rotation_angle = math.acos(max(-1.0, min(1.0, dot_product))) * 180 / math.pi
        
        if rot_length > 1e-6:
            transform.RotateWXYZ(rotation_angle, rot_x, rot_y, rot_z)
        
        transform.Scale(1, bond_length, 1)
        
        # Apply transform
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetInputConnection(cylinder.GetOutputPort())
        transform_filter.SetTransform(transform)
        
        # Create mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(transform_filter.GetOutputPort())
        
        # Set up actor
        self.actor.SetMapper(mapper)
        self.set_color(0.5, 0.5, 0.5)  # Gray bonds by default


class MoleculeRenderer:
    def __init__(self, renderer):
        self.renderer = renderer
        self.atoms = []
        self.bonds = []
        self.axes = None
        self.legend = None
    
    def clear(self):
        self.renderer.RemoveAllViewProps()
        self.atoms = []
        self.bonds = []
        self.axes = None
        self.legend = None
    
    def render_molecule(self, molecule_data):
        self.clear()
        
        atoms_data = molecule_data.get('atoms', [])
        bonds_data = molecule_data.get('bonds', [])
        
        if not atoms_data:
            return False
        
        # Render atoms
        for atom_info in atoms_data:
            atom = Atom(
                element=atom_info['element'],
                position=atom_info['position'],
                radius=0.3
            )
            self.atoms.append(atom)
            self.renderer.AddActor(atom.get_actor())
        
        # Render bonds
        for bond_info in bonds_data:
            bond = Bond(
                atom1_pos=bond_info['atom1_pos'],
                atom2_pos=bond_info['atom2_pos'],
                radius=0.1
            )
            self.bonds.append(bond)
            self.renderer.AddActor(bond.get_actor())
        
        # Add coordinate axes
        self.add_axes()
        
        # Add element legend/key
        self.add_legend(molecule_data)
        
        # Reset camera to fit molecule
        self.renderer.ResetCamera()
        
        return True
    
    def add_axes(self):
        self.axes = vtk.vtkAxesActor()
        self.axes.SetTotalLength(2, 2, 2)
        self.axes.SetShaftTypeToCylinder()
        self.renderer.AddActor(self.axes)
    
    def remove_axes(self):
        if self.axes:
            self.renderer.RemoveActor(self.axes)
            self.axes = None
    
    def add_legend(self, molecule_data):
        atoms_data = molecule_data.get('atoms', [])
        
        # Get unique elements in the molecule
        elements_present = {}
        for atom in atoms_data:
            element = atom['element']
            if element not in elements_present:
                elements_present[element] = 0
            elements_present[element] += 1
        
        # Sort elements by count (most common first)
        sorted_elements = sorted(elements_present.items(), key=lambda x: x[1], reverse=True)
        
        # Create legend actor
        self.legend = vtk.vtkLegendBoxActor()
        self.legend.SetNumberOfEntries(len(sorted_elements))
        
        # Element colours (same as CPK coloring)
        colours = {
            'H': (1.0, 1.0, 1.0),    # White
            'C': (0.5, 0.5, 0.5),    # Gray
            'N': (0.0, 0.0, 1.0),    # Blue
            'O': (1.0, 0.0, 0.0),    # Red
            'F': (0.0, 1.0, 0.0),    # Green
            'Cl': (0.0, 1.0, 0.0),   # Green
            'Br': (0.6, 0.2, 0.0),   # Brown
            'I': (0.5, 0.0, 0.5),    # Purple
            'S': (1.0, 1.0, 0.0),    # Yellow
            'P': (1.0, 0.5, 0.0),    # Orange
        }
        
        # Add each element to legend
        for idx, (element, count) in enumerate(sorted_elements):
            # Create a small sphere for the icon
            sphere = vtk.vtkSphereSource()
            sphere.SetRadius(0.5)
            sphere.SetPhiResolution(10)
            sphere.SetThetaResolution(10)
            
            # Get colour for this element
            colour = colours.get(element, (1.0, 0.0, 1.0))  # Magenta for unknown
            
            # Set legend entry
            label = f"{element} ({count})"
            self.legend.SetEntry(idx, sphere.GetOutput(), label, colour)
        
        # Position and style the legend
        self.legend.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        self.legend.GetPositionCoordinate().SetValue(0.02, 0.02)  # Bottom-left corner
        self.legend.SetWidth(0.15)
        self.legend.SetHeight(0.08 * len(sorted_elements))
        
        # Style the legend box
        self.legend.UseBackgroundOn()
        self.legend.SetBackgroundColor(0.2, 0.2, 0.2)
        self.legend.SetBackgroundOpacity(0.8)
        self.legend.BorderOff()
        
        # Add to renderer
        self.renderer.AddActor(self.legend)
    
    def remove_legend(self):
        if self.legend:
            self.renderer.RemoveActor(self.legend)
            self.legend = None
    
    def get_atom_count(self):
        return len(self.atoms)
    
    def get_bond_count(self):
        return len(self.bonds)