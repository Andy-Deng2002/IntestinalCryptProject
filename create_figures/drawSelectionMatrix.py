import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon, ConnectionPatch
import numpy as np

def create_hex_grid_data(rows, cols, size=1.0):
    polygons = []
    centers = []
    
    # Hexagon geometry (Pointy topped)
    dx = size * np.sqrt(3)
    dy = size * 1.5
    
    # We want Row 2 at top, Row 0 at bottom visually
    # Let's plot r=0 at y=0 (Bottom). r=2 at y=2*dy (Top).
    
    for r in range(rows):
        for c in range(cols):
            # Calculate center
            x_off = (dx / 2) if (r % 2 == 1) else 0
            cx = c * dx + x_off
            cy = r * dy
            
            centers.append((cx, cy))
            
            # Vertices
            angles = np.radians([30, 90, 150, 210, 270, 330])
            
            hx = cx + size * np.cos(angles)
            hy = cy + size * np.sin(angles)
            
            polygons.append(np.column_stack((hx, hy)))
            
    return polygons, centers, dx, dy

# Set up the figure
fig = plt.figure(figsize=(16, 8))

# Layout: 
# Left: 2D Hex Grid (Centered vertically)
# Right Top: Flattened Vector
# Right Bottom: Adjacency Matrix

gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1.2], height_ratios=[1, 1.5], 
                      wspace=0.25, hspace=0.05)

# Left subplot spans both rows
ax_grid = fig.add_subplot(gs[:, 0])

# Right subplots
ax_flat = fig.add_subplot(gs[0, 1])
ax_adj = fig.add_subplot(gs[1, 1])

# Parameters
rows = 3
cols = 4
n_cells = rows * cols
hex_size = 0.6

# Data preparation
polygons, centers, dx, dy = create_hex_grid_data(rows, cols, size=hex_size)

# Define States (0: WT, 1: MT)
states = np.zeros(n_cells, dtype=int)
# Set some mutants
# Manual 1D indices mapping: idx = r * cols + c
states[[1, 6, 8, 11]] = 1 

colors = {0: '#E6E6FA', 1: '#FF6347'} # Lavender (WT), Tomato (MT)

# --- Draw 2D Hex Grid (Left) ---
ax_grid.set_aspect('equal')
ax_grid.axis('off')
ax_grid.set_title("3x4 Hexagonal Offset Grid", fontsize=14, y=0.9)

# Draw Hexagons
for i, poly_verts in enumerate(polygons):
    c_state = states[i]
    c_color = colors[c_state]
    
    # Polygon
    poly = Polygon(poly_verts, facecolor=c_color, edgecolor='black', linewidth=1.0)
    ax_grid.add_patch(poly)
    
    cx, cy = centers[i]
    
    # Index Text (center)
    ax_grid.text(cx, cy - 0.25, str(i), ha='center', va='center', fontsize=8)
    
    # Prolif Text (below center)
    label = f'MT (rate {r'$\lambda$'})' if c_state == 1 else f'WT (rate 1)'
    ax_grid.text(cx, cy + 0.1, label, ha='center', va='center', fontsize=8, fontweight='bold')

# Add Row Selection Policy Labels (Alpha)
# Top Row (Highest r) -> Alpha
# Bottom Row (Lowest r) -> Alpha^3
row_alphas_labels = [r'$\alpha^3$', r'$\alpha^2$', r'$\alpha$'] # r=0, r=1, r=2

for r in range(rows):
    # Position: Left of the row
    # Use r to find y
    cy = r * dy
    # x should be to the left of column 0
    # col 0 x is 0 (+ offset)
    x_label = -0.8
    
    label = row_alphas_labels[r]
    ax_grid.text(x_label, cy, label, ha='center', va='center', fontsize=14, fontweight='bold')

# Bounding box for grid
ax_grid.set_xlim(-2, cols * dx)
ax_grid.set_ylim(-dy, rows * dy + 0.5)


# --- Draw Flattened Vector (Top Right) ---
ax_flat.set_xlim(-1, n_cells)
ax_flat.set_ylim(-2.5, 1.5)
ax_flat.set_aspect('equal', adjustable='datalim')
ax_flat.axis('off')
ax_flat.set_title("Flattened Representation", fontsize=14, y=0.75)
cell_radius = 0.4
y_center = -0.3  # Shift down

for i in range(n_cells):
    c_color = colors[states[i]]
    # Draw circle
    circle = patches.Circle((i, y_center), radius=cell_radius, edgecolor='black', facecolor=c_color, linewidth=1.5)
    ax_flat.add_patch(circle)
    # Label inside
    ax_flat.text(i, y_center, str(states[i]), ha='center', va='center', fontsize=10, fontweight='bold')
    # Index above
    ax_flat.text(i, y_center + 0.8, str(i), ha='center', va='center', fontsize=10, color='gray')

# Draw Pi Vector below
# Map indices to alphas
# r=2 (11, 10, 9, 8) -> Alpha ? No.
# Indices 8,9,10,11 are top. They get Alpha.
# Indices 0,1,2,3 are bottom. They get Alpha^3.

pi_mapping = {0: r'$\alpha^3$', 1: r'$\alpha^2$', 2: r'$\alpha$'}

ax_flat.text(-0.5, y_center - 1.2, r'$\vec{\pi} = ($', ha='right', va='center', fontsize=14)
for i in range(n_cells):
    r_idx = i // cols
    label = pi_mapping[r_idx]
    ax_flat.text(i, y_center - 1.2, label, ha='center', va='center', fontsize=11)
    if i < n_cells - 1:
        ax_flat.text(i + 0.5, y_center - 1.2, ",", ha='center', va='center', fontsize=12)
ax_flat.text(n_cells - 0.6, y_center - 1.2, r'$)$', ha='left', va='center', fontsize=14)


# --- Draw Adjacency Matrix (Bottom Right) ---
ax_adj.set_title("Adjacency Matrix", fontsize=14)

adj_matrix = np.zeros((n_cells, n_cells), dtype=int)

# Neighbor logic
# "Periodic in X, Discrete/Fixed in Y"
directions_even = [
    (0, -1), (0, 1),
    (-1, -1), (-1, 0),
    (1, -1), (1, 0)
]
directions_odd = [
    (0, -1), (0, 1),
    (-1, 0), (-1, 1),
    (1, 0), (1, 1)
]

for r in range(rows):
    for c in range(cols):
        u_idx = r * cols + c
        
        # Directions depend on parity of r
        dirs = directions_odd if (r % 2 == 1) else directions_even
        
        for dr, dc in dirs:
            nr, nc = r + dr, c + dc
            
            # Check Y bounds (Fixed)
            if 0 <= nr < rows:
                # Check X bounds (Periodic)
                nc = nc % cols
                v_idx = nr * cols + nc
                adj_matrix[u_idx, v_idx] = 1

ax_adj.imshow(adj_matrix, cmap='Purples', interpolation='nearest', aspect='equal')
# Ticks
ax_adj.set_xticks(np.arange(n_cells))
ax_adj.set_yticks(np.arange(n_cells))
ax_adj.set_xticklabels(np.arange(n_cells), fontsize=8)
ax_adj.set_yticklabels(np.arange(n_cells), fontsize=8)

# Grid
ax_adj.set_xticks(np.arange(n_cells + 1) - 0.5, minor=True)
ax_adj.set_yticks(np.arange(n_cells + 1) - 0.5, minor=True)
ax_adj.grid(which="minor", color="black", linestyle='-', linewidth=1)
ax_adj.tick_params(which="minor", size=0)
ax_adj.tick_params(which="major", labelsize=8)


# --- Connection Arrows ---

# 1. Grid (Left) -> Vector (Top Right)
# Use figure coordinates for better control
from matplotlib.patches import FancyArrowPatch

# Arrow 1: From right side of grid to left side of vector
arrow1 = FancyArrowPatch((0.45, 0.64), (0.55, 0.73),
                        transform=fig.transFigure,
                        arrowstyle='-|>', mutation_scale=20, 
                        linewidth=1.5, linestyle='--',
                        color='black', zorder=10)
fig.patches.append(arrow1)

# Arrow 2: From right side of grid to left side of matrix
arrow2 = FancyArrowPatch((0.47, 0.35), (0.58, 0.35),
                        transform=fig.transFigure,
                        arrowstyle='-|>', mutation_scale=20, 
                        linewidth=1.5, linestyle='--',
                        color='black', zorder=10)
fig.patches.append(arrow2)

plt.tight_layout()
plt.savefig('micSMM_illustration.png', dpi=300, bbox_inches='tight', transparent=True)
plt.show()