import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# Grid settings
rows = 4
cols = 6
cell_size = 1.0

# Define cell types: 0 for Wild-type, 1 for Mutant
# Let's create a random looking distribution or a specific pattern
grid_data = np.zeros((rows, cols), dtype=int)

# Set some specific cells as mutants (creating a cluster typically seen in expansion)
mutant_locations = [(1, 2), (1, 3), (2, 2), (2, 3), (2,1)]
for r, c in mutant_locations:
    grid_data[r, c] = 1

# Colors
color_wt = '#e6e6e6'  # Light grey for Wild-type
color_mut = '#FF6B6B' # Soft red for Mutant
text_color = 'black'

# Create figure
fig, ax = plt.subplots(figsize=(6, 4))

# Draw the grid
for r in range(rows):
    for c in range(cols):
        # Determine color and label based on type
        if grid_data[r, c] == 1:
            face_color = color_mut
            label = r'$r$' # Proliferative advantage r or lambda
            
            # Optional: Add a subtle glow or darker border for mutants
            edge_color = '#c92a2a'
            line_width = 2
        else:
            face_color = color_wt
            label = '$1$'
            edge_color = 'white'
            line_width = 1.5
        
        # Create a rectangle (cell)
        # Note: matplotlib coords (x,y) correspond to (col, row)
        # We invert row index so row 0 is at top
        rect = patches.Rectangle(
            (c * cell_size, (rows - 1 - r) * cell_size), 
            cell_size, cell_size, 
            linewidth=line_width, 
            edgecolor=edge_color, 
            facecolor=face_color
        )
        ax.add_patch(rect)
        
        # Add text label in center
        cx = c * cell_size + cell_size / 2
        cy = (rows - 1 - r) * cell_size + cell_size / 2
        ax.text(cx, cy, label, ha='center', va='center', fontsize=20, color=text_color, fontweight='bold')

# Set limits and aspect
ax.set_xlim(0, cols * cell_size)
ax.set_ylim(0, rows * cell_size)
ax.set_aspect('equal')

# Remove axes ticks and labels for a clean look
ax.axis('off')

# Add a title or legend-like annotation if desired (Optional)
# plt.title("Stochastic Spatial Model: Grid Topology", fontsize=14)

# Save the figure with transparent background
plt.tight_layout()
plt.savefig('cell_grid.png', dpi=300, transparent=True)
plt.show()