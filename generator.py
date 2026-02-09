import numpy as np
import math


# Atomic data (atomic radius in angstroms, color)
ATOM_DATA = {
    'H': {'radius': 0.31, 'color': (255, 255, 255), 'name': 'Hydrogen'},
    'C': {'radius': 0.76, 'color': (50, 50, 50), 'name': 'Carbon'},
    'N': {'radius': 0.71, 'color': (50, 100, 255), 'name': 'Nitrogen'},
    'O': {'radius': 0.66, 'color': (255, 50, 50), 'name': 'Oxygen'},
    'F': {'radius': 0.57, 'color': (144, 224, 144), 'name': 'Fluorine'},
    'P': {'radius': 1.07, 'color': (255, 128, 0), 'name': 'Phosphorus'},
    'S': {'radius': 1.05, 'color': (255, 255, 0), 'name': 'Sulfur'},
    'Cl': {'radius': 1.02, 'color': (50, 255, 50), 'name': 'Chlorine'},
    'Br': {'radius': 1.20, 'color': (165, 42, 42), 'name': 'Bromine'},
    'I': {'radius': 1.39, 'color': (148, 0, 211), 'name': 'Iodine'},
}


class Matrix3D:
    @staticmethod
    def rotation_x(angle):
        # Rotation matrix around X axis
        c, s = np.cos(angle), np.sin(angle)
        return np.array([
            [1, 0, 0],
            [0, c, -s],
            [0, s, c]
        ])
    
    @staticmethod
    def rotation_y(angle):
        # Rotation matrix around Y axis
        c, s = np.cos(angle), np.sin(angle)
        return np.array([
            [c, 0, s],
            [0, 1, 0],
            [-s, 0, c]
        ])
    
    @staticmethod
    def rotation_z(angle):
        # Rotation matrix around Z axis
        c, s = np.cos(angle), np.sin(angle)
        return np.array([
            [c, -s, 0],
            [s, c, 0],
            [0, 0, 1]
        ])
    
    @staticmethod
    def project_3d_to_2d(point, fov, viewer_distance, window_width, window_height, ui_panel_width):
        x, y, z = point
        factor = fov / (viewer_distance + z)
        x_proj = x * factor + window_width // 2 - ui_panel_width // 2
        y_proj = -y * factor + window_height // 2
        return int(x_proj), int(y_proj), z


class Atom:
    def __init__(self, element, position, atom_id):
        self.element = element
        self.position = np.array(position, dtype=float)
        self.atom_id = atom_id
        self.data = ATOM_DATA.get(element, ATOM_DATA['C'])
    
    def get_screen_pos(self, rotation_matrix, fov, viewer_distance, scale, window_width, window_height, ui_panel_width):
        rotated = rotation_matrix @ (self.position * scale)
        return Matrix3D.project_3d_to_2d(rotated, fov, viewer_distance, window_width, window_height, ui_panel_width)


class Bond:
    # Represents a bond between two atoms
    
    def __init__(self, atom1, atom2, bond_order=1):
        self.atom1 = atom1
        self.atom2 = atom2
        self.bond_order = bond_order  # 1=single, 2=double, 3=triple


class Molecule:
    def __init__(self, name, iupac_name=None):
        self.name = name
        self.iupac_name = iupac_name or name
        self.atoms = []
        self.bonds = []
    
    def add_atom(self, element, position):
        atom_id = len(self.atoms) + 1
        atom = Atom(element, position, atom_id)
        self.atoms.append(atom)
        return atom
    
    def add_bond(self, atom1, atom2, bond_order=1):
        bond = Bond(atom1, atom2, bond_order)
        self.bonds.append(bond)
    
    def calculate_distances(self):
        # Calculate distance matrix between all atoms
        n = len(self.atoms)
        distance_matrix = np.zeros((n, n))
        
        for i, atom1 in enumerate(self.atoms):
            for j, atom2 in enumerate(self.atoms):
                if i != j:
                    dist = np.linalg.norm(atom1.position - atom2.position)
                    distance_matrix[i][j] = dist
        
        return distance_matrix
    
    def to_dict(self):
        return {
            'name': self.name,
            'iupac_name': self.iupac_name,
            'atoms': [
                {
                    'id': atom.atom_id,
                    'element': atom.element,
                    'position': atom.position.tolist()
                }
                for atom in self.atoms
            ],
            'bonds': [
                {
                    'atom1_id': bond.atom1.atom_id,
                    'atom2_id': bond.atom2.atom_id,
                    'order': bond.bond_order
                }
                for bond in self.bonds
            ]
        }
    
    @classmethod
    def from_dict(cls, data):
        mol = cls(data['name'], data.get('iupac_name'))
        
        # Create atoms
        atom_map = {}
        for atom_data in data['atoms']:
            atom = mol.add_atom(
                atom_data['element'],
                atom_data['position']
            )
            atom_map[atom_data['id']] = atom
        
        # Create bonds
        for bond_data in data['bonds']:
            atom1 = atom_map[bond_data['atom1_id']]
            atom2 = atom_map[bond_data['atom2_id']]
            mol.add_bond(atom1, atom2, bond_data.get('order', 1))
        
        return mol


class MoleculeLibrary:
    # Predefined molecules with accurate 3D coordinates
    @staticmethod
    def create_methane():
        mol = Molecule("Methane", "methane")
        
        # Central carbon
        c = mol.add_atom('C', [0, 0, 0])
        
        # Tetrahedral hydrogen positions
        angle = math.radians(109.47)  # Tetrahedral angle
        bond_length = 1.09
        
        positions = [
            [bond_length, 0, 0],
            [-bond_length * math.cos(angle), bond_length * math.sin(angle), 0],
            [-bond_length * math.cos(angle), -bond_length * math.sin(angle) * math.cos(math.radians(60)), 
             bond_length * math.sin(angle) * math.sin(math.radians(60))],
            [-bond_length * math.cos(angle), -bond_length * math.sin(angle) * math.cos(math.radians(60)), 
             -bond_length * math.sin(angle) * math.sin(math.radians(60))],
        ]
        
        for pos in positions:
            h = mol.add_atom('H', pos)
            mol.add_bond(c, h)
        
        return mol
    
    @staticmethod
    def create_ethanol():
        mol = Molecule("Ethanol", "ethanol")
        
        # C-C backbone
        c1 = mol.add_atom('C', [-0.77, 0, 0])
        c2 = mol.add_atom('C', [0.77, 0, 0])
        mol.add_bond(c1, c2)
        
        # OH group
        o = mol.add_atom('O', [1.2, 1.2, 0])
        h_oh = mol.add_atom('H', [2.0, 1.4, 0])
        mol.add_bond(c2, o)
        mol.add_bond(o, h_oh)
        
        # Hydrogens on C1
        h1 = mol.add_atom('H', [-1.2, 0.8, 0.8])
        h2 = mol.add_atom('H', [-1.2, -0.8, 0.8])
        h3 = mol.add_atom('H', [-1.2, 0, -1.2])
        mol.add_bond(c1, h1)
        mol.add_bond(c1, h2)
        mol.add_bond(c1, h3)
        
        # Hydrogens on C2
        h4 = mol.add_atom('H', [1.2, -0.6, 0.8])
        h5 = mol.add_atom('H', [1.2, -0.6, -0.8])
        mol.add_bond(c2, h4)
        mol.add_bond(c2, h5)
        
        return mol
    
    @staticmethod
    def create_water():
        mol = Molecule("Water", "oxidane")
        
        o = mol.add_atom('O', [0, 0, 0])
        
        # Bent geometry (104.5 degrees)
        angle = math.radians(104.5)
        bond_length = 0.96
        
        h1 = mol.add_atom('H', [bond_length * math.sin(angle/2), bond_length * math.cos(angle/2), 0])
        h2 = mol.add_atom('H', [-bond_length * math.sin(angle/2), bond_length * math.cos(angle/2), 0])
        
        mol.add_bond(o, h1)
        mol.add_bond(o, h2)
        
        return mol
    
    @staticmethod
    def create_benzene():
        mol = Molecule("Benzene", "benzene")
        
        # Hexagonal ring
        radius = 1.40  # C-C bond length in benzene
        carbons = []
        
        for i in range(6):
            angle = math.radians(60 * i)
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            c = mol.add_atom('C', [x, y, 0])
            carbons.append(c)
        
        # Carbon-carbon bonds (alternating single/double)
        for i in range(6):
            bond_order = 2 if i % 2 == 0 else 1
            mol.add_bond(carbons[i], carbons[(i + 1) % 6], bond_order)
        
        # Add hydrogens
        h_dist = 1.09
        for i, c in enumerate(carbons):
            angle = math.radians(60 * i)
            x = (radius + h_dist) * math.cos(angle)
            y = (radius + h_dist) * math.sin(angle)
            h = mol.add_atom('H', [x, y, 0])
            mol.add_bond(c, h)
        
        return mol
    
    @staticmethod
    def create_ammonia():
        mol = Molecule("Ammonia", "azane")
        
        n = mol.add_atom('N', [0, 0, 0])
        
        # Trigonal pyramidal geometry
        bond_length = 1.01
        angle = math.radians(107)
        
        h1 = mol.add_atom('H', [bond_length, 0, 0])
        h2 = mol.add_atom('H', [-bond_length * math.cos(angle), bond_length * math.sin(angle), 0])
        h3 = mol.add_atom('H', [-bond_length * math.cos(angle), -bond_length * math.sin(angle) * 0.5, 
                                 bond_length * math.sin(angle) * 0.866])
        
        mol.add_bond(n, h1)
        mol.add_bond(n, h2)
        mol.add_bond(n, h3)
        
        return mol