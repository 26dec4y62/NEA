from rdkit import Chem

class Molecule:
    def __init__(self):
        pass
    ATOM_COLOURS = {
        "H": (1.0, 1.0, 1.0),
        "C": (0.2, 0.2, 0.2),
        "N": (0.0, 0.0, 1.0),
        "O": (1.0, 0.0, 0.0),
        "F": (0.0, 0.8, 0.0),
        "Cl": (0.0, 1.0, 0.0),
        "Br": (0.6, 0.2, 0.2),
        "I": (0.4, 0.0, 0.6),
        "S": (1.0, 0.8, 0.2),
        "P": (1.0, 0.5, 0.0)
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

    @classmethod
    def structural_formula(cls, smiles):
        # Heuristic algorithm
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")

        Chem.Kekulize(mol, clearAromaticFlags=True)
        visited = set()
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()

        # Group fused rings into systems
        ring_systems = []
        used = set()
        for i, ring in enumerate(atom_rings):
            if i in used:
                continue
            fused = set(ring)
            for j, other in enumerate(atom_rings):
                if i == j or j in used:
                    continue
                if fused.intersection(other):
                    fused.update(other)
                    used.add(j)
            used.add(i)
            ring_systems.append(fused)
        atom_to_system = {}
        for si, sys in enumerate(ring_systems):
            for a in sys:
                atom_to_system[a] = si

        def atom_to_group(atom):
            sym = atom.GetSymbol()
            hcount = atom.GetTotalNumHs()
            if sym == "C":
                return f"C{'H'+str(hcount) if hcount > 0 else ''}"
            elif sym == "O":
                return "OH" if hcount == 1 else "O"
            elif sym == "N":
                return "NH2" if hcount == 2 else "N"
            else:
                return sym

        def handle_ring_system(system_atoms):
            atoms = [mol.GetAtomWithIdx(i) for i in system_atoms]
            carbons = [a for a in atoms if a.GetSymbol() == "C"]
            hydrogens = sum(a.GetTotalNumHs() for a in atoms)
            hetero = [a for a in atoms if a.GetSymbol() not in ("C", "H")]

            if carbons and all(a.GetIsAromatic() for a in carbons) and not hetero:
                return f"C{len(carbons)}H{hydrogens}"
            return "".join(atom_to_group(a) for a in atoms)

        def walk(atom, prev=None):
            idx = atom.GetIdx()
            if idx in visited:
                return ""
            visited.add(idx)

            # ring system handling
            if idx in atom_to_system:
                sys_idx = atom_to_system[idx]
                sys_atoms = ring_systems[sys_idx]
                if not visited.intersection(sys_atoms):
                    visited.update(sys_atoms)
                    return handle_ring_system(sys_atoms)

            group = atom_to_group(atom)
            branches = []

            for nbr in atom.GetNeighbors():
                if prev is not None and nbr.GetIdx() == prev.GetIdx():
                    continue
                if nbr.GetIdx() in visited:
                    continue
                branch_str = walk(nbr, atom)
                if branch_str:
                    branches.append(branch_str)

            if len(branches) == 0:
                return group
            elif len(branches) == 1:
                return group + branches[0]
            else:
                main = branches[0]
                side_branches = "".join([f"({b})" for b in branches[1:]])
                return group + main + side_branches

        # start at a carbon if possible
        start = next((a for a in mol.GetAtoms() if a.GetSymbol() == "C"), mol.GetAtomWithIdx(0))
        return walk(start)

    def _extract_atoms_and_bonds_rdkit(self):
        if not self.self.mol:
            return
        
        # Get conformer (3D coordinates)
        conf = self.self.mol.GetConformer()
        
        # Extract atoms
        self.atoms = []
        for atom in self.self.mol.GetAtoms():
            idx = atom.GetIdx()
            pos = conf.GetAtomPosition(idx)
            self.atoms.append({
                'symbol': atom.GetSymbol(),
                'position': [pos.x, pos.y, pos.z],
                'index': idx
            })
        
        # Extract bonds
        self.bonds = []
        for bond in self.self.mol.GetBonds():
            self.bonds.append({
                'atoms': [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()],
                'order': bond.GetBondType()
            })
    