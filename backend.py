import json
from urllib.parse import quote
from urllib.request import urlopen
from urllib.error import URLError, HTTPError
from collections import deque
from collections import Counter

FILE = "molecules.json"

class appdata:
    def __init__(self):
        self.history = []
        self.molecule_cache = {}  # Cache molecules by name
        self.load_data()

    def load_data(self):
        try:
            with open(FILE, "r") as f:
                data = json.load(f)
            self.history = data.get("History", [])
            self.molecule_cache = data.get("Cache", {})
        except FileNotFoundError:
            self.history = []
            self.molecule_cache = {}
            self.save_data()

    def save_data(self):
        with open(FILE, "w") as f:
            json.dump({
                "History": self.history, 
                "Cache": self.molecule_cache
            }, f, indent=2)

    def add_history(self, name):
        if name not in self.history:
            self.history.append(name)
            self.save_data()

    def get_molecule(self, name):
        name_lower = name.lower()
        
        # Check cache first
        if name_lower in self.molecule_cache:
            return self.molecule_cache[name_lower]
        
        # Fetch from PubChem
        molecule_data = self.fetch_molecule(name)
        
        # Add to history and cache
        self.add_history(name)
        self.molecule_cache[name_lower] = molecule_data
        self.save_data()
        
        return molecule_data

    # Needs to be redone
    def fetch_molecule(self, name):
        # Fetch molecule data from PubChem including 3D structure
        name_enc = quote(name)
        
        # Get CID (Compound ID)
        cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name_enc}/cids/JSON"
        try:
            with urlopen(cid_url, timeout=10) as response:
                cid_data = json.loads(response.read().decode("utf-8"))
            cid = cid_data["IdentifierList"]["CID"][0]
        except (URLError, HTTPError, KeyError) as e:
            raise ValueError(f"Could not find molecule: {name}")
        
        # Get basic properties
        props_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName,MolecularFormula,IsomericSMILES/JSON"
        )
        try:
            with urlopen(props_url, timeout=10) as response:
                props_data = json.loads(response.read().decode("utf-8"))
            props = props_data["PropertyTable"]["Properties"][0]
        except (URLError, HTTPError, KeyError) as e:
            raise ValueError(f"Could not fetch properties for: {name}")
        
        # Get 3D structure - SDF format
        sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        try:
            with urlopen(sdf_url, timeout=10) as response:
                sdf_text = response.read().decode("utf-8")
            atoms, bonds = self.parse_sdf(sdf_text)
        except (URLError, HTTPError):
            # 2D coordinates fallback
            print(f"3D structure not available for {name}, using 2D")
            sdf_url_2d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
            try:
                with urlopen(sdf_url_2d, timeout=10) as response:
                    sdf_text = response.read().decode("utf-8")
                atoms, bonds = self.parse_sdf(sdf_text)
            except (URLError, HTTPError) as e:
                raise ValueError(f"Could not fetch structure for: {name}")
        
        # Generate structural formula
        structural_formula = self.generate_structural_formula(atoms, bonds)
        
        return {
            "cid": cid,
            "name": name,
            "iupac_name": props.get("IUPACName", "N/A"),
            "formula": props.get("MolecularFormula", "N/A"),
            "structural_formula": structural_formula,
            "smiles": props.get("IsomericSMILES", "N/A"),
            "atoms": atoms,
            "bonds": bonds
        }

    def parse_sdf(self, sdf_text):
        lines = sdf_text.split("\n")
        
        counts_line = lines[3].split()
        num_atoms = int(counts_line[0])
        num_bonds = int(counts_line[1])
        
        atoms = []
        for i in range(4, 4 + num_atoms):
            parts = lines[i].split()
            if len(parts) >= 4:
                atoms.append({
                    "element": parts[3],
                    "position": (float(parts[0]), float(parts[1]), float(parts[2]))
                })
        
        # Parse bonds
        bonds = []
        bond_start = 4 + num_atoms
        for i in range(bond_start, bond_start + num_bonds):
            parts = lines[i].split()
            if len(parts) >= 3:
                atom1_idx = int(parts[0]) - 1
                atom2_idx = int(parts[1]) - 1
                bond_type = int(parts[2])  # 1=single, 2=double, 3=triple
                
                bonds.append({
                    "atom1": atom1_idx,
                    "atom2": atom2_idx,
                    "atom1_pos": atoms[atom1_idx]["position"],
                    "atom2_pos": atoms[atom2_idx]["position"],
                    "type": bond_type
                })
        
        return atoms, bonds

    def generate_structural_formula(self, atoms, bonds):
        if not atoms:
            return "N/A"
        
        # Build adjacency list
        adjacency = {i: [] for i in range(len(atoms))}
        for bond in bonds:
            atom1 = bond["atom1"]
            atom2 = bond["atom2"]
            bond_type = bond["type"]
            adjacency[atom1].append((atom2, bond_type))
            adjacency[atom2].append((atom1, bond_type))
        
        # Find the "main chain" - typically the longest carbon chain
        main_chain = self._find_main_chain(atoms, adjacency)
        
        if not main_chain:
            # Fallback: just show molecular formula
            return self._generate_molecular_formula(atoms)
        
        # Build structural formula string
        formula_parts = []
        visited = set()
        
        for i, atom_idx in enumerate(main_chain):
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            
            element = atoms[atom_idx]["element"]
            
            # Get branches attached to this atom
            branches = []
            for neighbor_idx, bond_type in adjacency[atom_idx]:
                if neighbor_idx not in visited and neighbor_idx not in main_chain:
                    branch_formula = self._build_branch(neighbor_idx, atom_idx, atoms, adjacency, visited)
                    if branch_formula:
                        branches.append(branch_formula)
            
            # Add the main atom
            atom_str = element
            
            # Add hydrogens if it"s a carbon in the main chain
            if element == "C":
                h_count = self._count_hydrogens(atom_idx, atoms, adjacency, visited)
                if h_count > 0:
                    atom_str += f"H{h_count}" if h_count > 1 else "H"
            
            # Add branches in parentheses
            for branch in branches:
                atom_str += f"({branch})"
            
            formula_parts.append(atom_str)
            
            # Add bond symbol between atoms
            if i < len(main_chain) - 1:
                next_idx = main_chain[i + 1]
                bond_type = self._get_bond_type(atom_idx, next_idx, adjacency)
                if bond_type == 2:
                    formula_parts.append("=")
                elif bond_type == 3:
                    formula_parts.append("≡")
                # Single bonds are implied, no symbol needed
        
        result = "".join(formula_parts)
        
        # If result is too complex or empty, fall back to molecular formula
        if not result or len(result) > 50:
            return self._generate_molecular_formula(atoms)
        
        return result

    def _find_main_chain(self, atoms, adjacency):
        # Find the longest carbon chain in the molecule
        carbon_indices = [i for i, atom in enumerate(atoms) if atom["element"] == "C"]
        
        if not carbon_indices:
            # No carbons, just return first few atoms
            return list(range(min(3, len(atoms))))
        
        # Use BFS to find longest path through carbons
        longest_chain = []
        
        for start in carbon_indices:
            chain = self._bfs_longest_path(start, atoms, adjacency)
            if len(chain) > len(longest_chain):
                longest_chain = chain
        
        return longest_chain if longest_chain else carbon_indices[:1]

    def _bfs_longest_path(self, start, atoms, adjacency):
        queue = deque([(start, [start])])
        longest = [start]
        
        while queue:
            current, path = queue.popleft()
            
            if len(path) > len(longest):
                longest = path
            
            for neighbor, _ in adjacency[current]:
                if neighbor not in path and atoms[neighbor]["element"] == "C":
                    queue.append((neighbor, path + [neighbor]))
        return longest

    def _build_branch(self, atom_idx, parent_idx, atoms, adjacency, visited):
        # Recursive branch building
        if atom_idx in visited:
            return ""
        
        visited.add(atom_idx)
        element = atoms[atom_idx]["element"]
        
        # Count hydrogens on this atom
        h_count = self._count_hydrogens(atom_idx, atoms, adjacency, visited)
        
        branch = element
        if element == "C" and h_count > 0:
            branch += f"H{h_count}" if h_count > 1 else "H"
        
        # Recursively add sub-branches
        for neighbor_idx, bond_type in adjacency[atom_idx]:
            if neighbor_idx != parent_idx and neighbor_idx not in visited:
                sub_branch = self._build_branch(neighbor_idx, atom_idx, atoms, adjacency, visited)
                if sub_branch:
                    branch += f"({sub_branch})"
        
        return branch

    def _count_hydrogens(self, atom_idx, atoms, adjacency, visited):
        h_count = 0
        for neighbor_idx, _ in adjacency[atom_idx]:
            if neighbor_idx not in visited and atoms[neighbor_idx]["element"] == "H":
                h_count += 1
        return h_count

    def _get_bond_type(self, atom1_idx, atom2_idx, adjacency):
        for neighbor, bond_type in adjacency[atom1_idx]:
            if neighbor == atom2_idx:
                return bond_type
        return 1  # Default to single bond

    def _generate_molecular_formula(self, atoms):
        element_counts = Counter(atom["element"] for atom in atoms)
        formula_parts = []
        
        # Carbon first
        if "C" in element_counts:
            count = element_counts["C"]
            formula_parts.append(f"C{count}" if count > 1 else "C")
            del element_counts["C"]
        
        # Hydrogen second
        if "H" in element_counts:
            count = element_counts["H"]
            formula_parts.append(f"H{count}" if count > 1 else "H")
            del element_counts["H"]
        
        # Rest alphabetically
        for element in sorted(element_counts.keys()):
            count = element_counts[element]
            formula_parts.append(f"{element}{count}" if count > 1 else element)
        
        return "".join(formula_parts)

    def clear_cache(self):
        self.molecule_cache = {}
        self.save_data()

    def clear_history(self):
        self.history = []
        self.save_data()