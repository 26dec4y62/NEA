"""
Molecule Storage Module
Handles saving and loading molecules to/from JSON files
"""

import json
import os


class MoleculeStorage:
    """Backend storage operations for molecules"""
    
    def __init__(self, filename='saved_molecules.json'):
        self.filename = filename
        self.molecules = self._load_from_file()
    
    def _load_from_file(self):
        """Load previously saved molecules from JSON"""
        if os.path.exists(self.filename):
            try:
                with open(self.filename, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"Error loading molecules: {e}")
                return []
        return []
    
    def save_molecule(self, molecule_dict):
        """Save a molecule to the JSON file"""
        if not molecule_dict:
            return False
        
        # Check if already saved (update if exists)
        existing_index = None
        for i, m in enumerate(self.molecules):
            if m['name'] == molecule_dict['name']:
                existing_index = i
                break
        
        if existing_index is not None:
            # Update existing molecule
            self.molecules[existing_index] = molecule_dict
        else:
            # Add new molecule
            self.molecules.append(molecule_dict)
        
        # Write to file
        try:
            with open(self.filename, 'w') as f:
                json.dump(self.molecules, f, indent=2)
            return True
        except Exception as e:
            print(f"Error saving molecule: {e}")
            return False
    
    def get_saved_molecules(self):
        """Return list of all saved molecules"""
        return self.molecules.copy()
    
    def get_molecule_by_name(self, name):
        """Retrieve a specific molecule by name"""
        for mol in self.molecules:
            if mol['name'] == name:
                return mol
        return None
    
    def delete_molecule(self, name):
        """Delete a molecule by name"""
        self.molecules = [m for m in self.molecules if m['name'] != name]
        
        try:
            with open(self.filename, 'w') as f:
                json.dump(self.molecules, f, indent=2)
            return True
        except Exception as e:
            print(f"Error deleting molecule: {e}")
            return False