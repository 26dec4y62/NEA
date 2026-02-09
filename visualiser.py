import pygame
from pygame.locals import *
import math

from generator import (
    Molecule, ATOM_DATA, MoleculeLibrary, Matrix3D
)
from storage import MoleculeStorage


# Constants
WINDOW_WIDTH = 1200
WINDOW_HEIGHT = 800
FPS = 60
BACKGROUND_COLOR = (15, 15, 25)
UI_PANEL_WIDTH = 300


class MoleculeVisualizer:
    def __init__(self):
        pygame.init()
        self.screen = pygame.display.set_mode((WINDOW_WIDTH, WINDOW_HEIGHT))
        pygame.display.set_caption("3D Molecular Visualizer")
        self.clock = pygame.time.Clock()
        self.running = True
        
        # 3D viewing parameters
        self.rotation_x = 0
        self.rotation_y = 0
        self.rotation_z = 0
        self.scale = 80
        self.fov = 500
        self.viewer_distance = 5
        
        # Mouse control
        self.mouse_down = False
        self.last_mouse_pos = None
        
        # Storage backend
        self.storage = MoleculeStorage()
        
        # Current molecule
        self.molecule = None
        
        # UI state
        self.font = pygame.font.Font(None, 24)
        self.small_font = pygame.font.Font(None, 18)
        self.show_key = True
        
        # Load initial molecule
        self.load_molecule(MoleculeLibrary.create_methane())
    
    def load_molecule(self, molecule):
        self.molecule = molecule
        
        # Reset rotation
        self.rotation_x = 0
        self.rotation_y = 0
        self.rotation_z = 0
    
    def save_current_molecule(self):
        if self.molecule:
            mol_dict = self.molecule.to_dict()
            success = self.storage.save_molecule(mol_dict)
            if success:
                print(f"Saved {self.molecule.name}")
            return success
        return False
    
    def get_rotation_matrix(self):
        # Get combined rotation matrix
        rx = Matrix3D.rotation_x(self.rotation_x)
        ry = Matrix3D.rotation_y(self.rotation_y)
        rz = Matrix3D.rotation_z(self.rotation_z)
        return rz @ ry @ rx
    
    def draw_molecule(self):
        if not self.molecule:
            return
        
        rotation_matrix = self.get_rotation_matrix()
        
        # Get all atom screen positions with depth
        atom_positions = []
        for atom in self.molecule.atoms:
            x, y, z = atom.get_screen_pos(
                rotation_matrix, 
                self.fov, 
                self.viewer_distance, 
                self.scale,
                WINDOW_WIDTH,
                WINDOW_HEIGHT,
                UI_PANEL_WIDTH
            )
            atom_positions.append((atom, x, y, z))
        
        # Sort by depth (z-coordinate) for proper rendering
        atom_positions.sort(key=lambda p: p[3], reverse=True)
        
        # Draw bonds first
        for bond in self.molecule.bonds:
            atom1_data = next(p for p in atom_positions if p[0] == bond.atom1)
            atom2_data = next(p for p in atom_positions if p[0] == bond.atom2)
            
            x1, y1 = atom1_data[1], atom1_data[2]
            x2, y2 = atom2_data[1], atom2_data[2]
            
            # Draw bond lines
            if bond.bond_order == 1:
                pygame.draw.line(self.screen, (100, 100, 100), (x1, y1), (x2, y2), 3)
            elif bond.bond_order == 2:
                # Draw double bond as two parallel lines
                offset = 3
                dx = x2 - x1
                dy = y2 - y1
                length = math.sqrt(dx*dx + dy*dy)
                if length > 0:
                    perpx = -dy / length * offset
                    perpy = dx / length * offset
                    pygame.draw.line(self.screen, (100, 100, 100), 
                                   (x1 + perpx, y1 + perpy), (x2 + perpx, y2 + perpy), 2)
                    pygame.draw.line(self.screen, (100, 100, 100), 
                                   (x1 - perpx, y1 - perpy), (x2 - perpx, y2 - perpy), 2)
            elif bond.bond_order == 3:
                # Draw triple bond
                offset = 4
                dx = x2 - x1
                dy = y2 - y1
                length = math.sqrt(dx*dx + dy*dy)
                if length > 0:
                    perpx = -dy / length * offset
                    perpy = dx / length * offset
                    pygame.draw.line(self.screen, (100, 100, 100), (x1, y1), (x2, y2), 2)
                    pygame.draw.line(self.screen, (100, 100, 100), 
                                   (x1 + perpx, y1 + perpy), (x2 + perpx, y2 + perpy), 2)
                    pygame.draw.line(self.screen, (100, 100, 100), 
                                   (x1 - perpx, y1 - perpy), (x2 - perpx, y2 - perpy), 2)
        
        # Draw atoms
        for atom, x, y, z in atom_positions:
            # Scale radius based on depth for perspective
            depth_scale = self.fov / (self.viewer_distance + z)
            radius = int(atom.data['radius'] * 30 * depth_scale)
            radius = max(radius, 5)
            
            # Draw atom sphere
            pygame.draw.circle(self.screen, atom.data['color'], (x, y), radius)
            pygame.draw.circle(self.screen, (200, 200, 200), (x, y), radius, 1)
            
            # Draw atom ID
            id_text = self.small_font.render(str(atom.atom_id), True, (255, 255, 255))
            self.screen.blit(id_text, (x - id_text.get_width()//2, y - id_text.get_height()//2))
    
    def draw_ui(self):
        # UI Panel background
        panel_rect = pygame.Rect(WINDOW_WIDTH - UI_PANEL_WIDTH, 0, UI_PANEL_WIDTH, WINDOW_HEIGHT)
        pygame.draw.rect(self.screen, (30, 30, 40), panel_rect)
        pygame.draw.line(self.screen, (60, 60, 70), 
                        (WINDOW_WIDTH - UI_PANEL_WIDTH, 0), 
                        (WINDOW_WIDTH - UI_PANEL_WIDTH, WINDOW_HEIGHT), 2)
        
        y_offset = 20
        
        # Title
        if self.molecule:
            title = self.font.render("Molecule Viewer", True, (255, 255, 255))
            self.screen.blit(title, (WINDOW_WIDTH - UI_PANEL_WIDTH + 20, y_offset))
            y_offset += 40
            
            # Molecule name
            name_text = self.small_font.render(f"Name: {self.molecule.name}", True, (200, 200, 200))
            self.screen.blit(name_text, (WINDOW_WIDTH - UI_PANEL_WIDTH + 20, y_offset))
            y_offset += 25
            
            # IUPAC name
            iupac_text = self.small_font.render(f"IUPAC: {self.molecule.iupac_name}", True, (200, 200, 200))
            self.screen.blit(iupac_text, (WINDOW_WIDTH - UI_PANEL_WIDTH + 20, y_offset))
            y_offset += 35
            
            # Atom key
            if self.show_key:
                key_title = self.font.render("Atom Key:", True, (255, 255, 255))
                self.screen.blit(key_title, (WINDOW_WIDTH - UI_PANEL_WIDTH + 20, y_offset))
                y_offset += 30
                
                # Group atoms by element
                element_counts = {}
                for atom in self.molecule.atoms:
                    if atom.element not in element_counts:
                        element_counts[atom.element] = []
                    element_counts[atom.element].append(atom.atom_id)
                
                for element, ids in sorted(element_counts.items()):
                    # Draw color circle
                    color = ATOM_DATA[element]['color']
                    pygame.draw.circle(self.screen, color, 
                                     (WINDOW_WIDTH - UI_PANEL_WIDTH + 35, y_offset + 8), 8)
                    pygame.draw.circle(self.screen, (200, 200, 200), 
                                     (WINDOW_WIDTH - UI_PANEL_WIDTH + 35, y_offset + 8), 8, 1)
                    
                    # Element name and IDs
                    text = f"{ATOM_DATA[element]['name']} ({element})"
                    elem_text = self.small_font.render(text, True, (200, 200, 200))
                    self.screen.blit(elem_text, (WINDOW_WIDTH - UI_PANEL_WIDTH + 50, y_offset))
                    y_offset += 20
                    
                    # Show atom IDs
                    ids_str = f"  IDs: {', '.join(map(str, ids))}"
                    ids_text = self.small_font.render(ids_str, True, (150, 150, 150))
                    self.screen.blit(ids_text, (WINDOW_WIDTH - UI_PANEL_WIDTH + 50, y_offset))
                    y_offset += 25
        
        # Controls
        y_offset = WINDOW_HEIGHT - 200
        controls_title = self.font.render("Controls:", True, (255, 255, 255))
        self.screen.blit(controls_title, (WINDOW_WIDTH - UI_PANEL_WIDTH + 20, y_offset))
        y_offset += 30
        
        controls = [
            "Drag: Rotate molecule",
            "Scroll: Zoom in/out",
            "1-5: Load molecules",
            "S: Save molecule",
            "K: Toggle key",
            "ESC: Quit"
        ]
        
        for control in controls:
            text = self.small_font.render(control, True, (180, 180, 180))
            self.screen.blit(text, (WINDOW_WIDTH - UI_PANEL_WIDTH + 20, y_offset))
            y_offset += 20
    
    def handle_events(self):
        for event in pygame.event.get():
            if event.type == QUIT:
                self.running = False
            
            elif event.type == KEYDOWN:
                if event.key == K_ESCAPE:
                    self.running = False
                elif event.key == K_k:
                    self.show_key = not self.show_key
                elif event.key == K_s:
                    self.save_current_molecule()
                elif event.key == K_1:
                    self.load_molecule(MoleculeLibrary.create_methane())
                elif event.key == K_2:
                    self.load_molecule(MoleculeLibrary.create_water())
                elif event.key == K_3:
                    self.load_molecule(MoleculeLibrary.create_ammonia())
                elif event.key == K_4:
                    self.load_molecule(MoleculeLibrary.create_ethanol())
                elif event.key == K_5:
                    self.load_molecule(MoleculeLibrary.create_benzene())
            
            elif event.type == MOUSEBUTTONDOWN:
                if event.button == 1:  # Left click
                    if event.pos[0] < WINDOW_WIDTH - UI_PANEL_WIDTH:
                        self.mouse_down = True
                        self.last_mouse_pos = event.pos
                elif event.button == 4:  # Scroll up
                    self.scale = min(self.scale * 1.1, 200)
                elif event.button == 5:  # Scroll down
                    self.scale = max(self.scale * 0.9, 20)
            
            elif event.type == MOUSEBUTTONUP:
                if event.button == 1:
                    self.mouse_down = False
            
            elif event.type == MOUSEMOTION:
                if self.mouse_down and self.last_mouse_pos:
                    dx = event.pos[0] - self.last_mouse_pos[0]
                    dy = event.pos[1] - self.last_mouse_pos[1]
                    
                    self.rotation_y += dx * 0.01
                    self.rotation_x += dy * 0.01
                    
                    self.last_mouse_pos = event.pos
    
    def run(self):
        while self.running:
            self.handle_events()
            
            # Clear screen
            self.screen.fill(BACKGROUND_COLOR)
            
            # Draw molecule
            self.draw_molecule()
            
            # Draw UI
            self.draw_ui()
            
            # Update display
            pygame.display.flip()
            self.clock.tick(FPS)
        
        pygame.quit()


def main():
    app = MoleculeVisualizer()
    app.run()


if __name__ == "__main__":
    main()