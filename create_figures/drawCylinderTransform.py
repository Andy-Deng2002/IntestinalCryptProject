import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon, FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def create_hex_grid(rows, cols, size=1.0):
    polygons = []
    
    # Hexagon geometry (Pointy topped)
    dx = size * np.sqrt(3)
    dy = size * 1.5
    
    width = cols * dx
    height = (rows * 1.5 + 0.5) * size
    
    # Generate Hexagons
    for r in range(rows):
        for c in range(cols):
            x_off = (dx / 2) if (r % 2 == 1) else 0
            cx = c * dx + x_off
            cy = r * dy
            
            # Vertices
            angles = np.radians([30, 90, 150, 210, 270, 330])
            hx = cx + size * np.cos(angles)
            hy = cy + size * np.sin(angles)
            
            polygons.append(np.column_stack((hx, hy)))
            
    # Generate Zig-Zag Cut Line (Left boundary)
    # The boundary follows the left edge of the hexagonal grid
    # For pointy-top hexagons with offset rows, we trace:
    # - Even rows: left vertices at index 2 (top-left) and 3 (bottom-left)
    # - Odd rows: they share vertices with even rows
    
    cut_line_points = []
    
    for r in range(rows):
        c = 0
        x_off = (dx / 2) if (r % 2 == 1) else 0
        cx = c * dx + x_off
        cy = r * dy
        
        angles = np.radians([30, 90, 150, 210, 270, 330])
        hx = cx + size * np.cos(angles)
        hy = cy + size * np.sin(angles)
        
        # Trace the leftmost boundary
        # Even row (r%2==0): protrudes to the left
        # Vertices: 2 (150°, top-left), 3 (210°, bottom-left)
        # Odd row (r%2==1): recessed
        # Vertices: 2 (150°, top-left), 3 (210°, bottom-left)
        
        if r == 0:
            # Start from bottom-left of first hexagon
            cut_line_points.append([hx[3], hy[3]])
        
        # Add the top-left vertex
        cut_line_points.append([hx[2], hy[2]])
        
        # Connect to next row's bottom-left if needed
        if r < rows - 1:
            # Calculate next row's hexagon position
            next_r = r + 1
            next_x_off = (dx / 2) if (next_r % 2 == 1) else 0
            next_cx = 0 * dx + next_x_off
            next_cy = next_r * dy
            
            next_hx = next_cx + size * np.cos(angles)
            next_hy = next_cy + size * np.sin(angles)
            
            # Add connection to next row's top-left vertex
            cut_line_points.append([next_hx[3], next_hy[3]])
    
    # Convert points to line segments
    cut_line_segments = []
    for i in range(len(cut_line_points) - 1):
        cut_line_segments.append((np.array(cut_line_points[i]), np.array(cut_line_points[i+1])))
            
    return polygons, width, height, cut_line_segments

def wrap_point(x, y, width, radius, rotation_offset_theta=-np.pi/2):
    theta = (x / width) * 2 * np.pi
    theta += rotation_offset_theta
    X = radius * np.cos(theta)
    Y = radius * np.sin(theta)
    Z = y
    return X, Y, Z

def draw_figure():
    rows = 10
    cols = 16
    size = 0.5
    
    polygons, width, height, cut_segments = create_hex_grid(rows, cols, size)
    radius = width / (2 * np.pi)
    
    fig = plt.figure(figsize=(11, 4.5))
    
    # --- Left: 3D Cylinder ---
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.set_axis_off()
    ax1.set_position([0.01, 0.05, 0.45, 0.9])  # [left, bottom, width, height] in figure coordinates
    
    # Center the cut to the front
    # The cut is roughly at x ~ 0.
    # x=0 -> theta = 0. We want theta to be -pi/2 (front view usually).
    offset = -np.pi/2
    
    # Draw Cylinder Polygons
    # Use Poly3DCollection for occlusion
    verts_3d = []
    for p in polygons:
        xs, ys, zs = wrap_point(p[:,0], p[:,1], width, radius, offset)
        verts_3d.append(list(zip(xs, ys, zs)))
        
    # Facecolor needs to be opaque to hide back lines
    poly_col = Poly3DCollection(verts_3d, facecolor='white', edgecolor='#555555', linewidths=0.5, alpha=1.0)
    ax1.add_collection3d(poly_col)
    
    # Draw Zig-Zag Cut Line (Red) on 3D
    for p_start, p_end in cut_segments:
        # Wrap start
        xs1, ys1, zs1 = wrap_point(p_start[0], p_start[1], width, radius, offset)
        # Wrap end
        xs2, ys2, zs2 = wrap_point(p_end[0], p_end[1], width, radius, offset)
        
        # To avoid "cutting through" the cylinder for long segments (shouldn't happen for short hex edges), just plot line
        ax1.plot([xs1, xs2], [ys1, ys2], [zs1, zs2], color='red', linewidth=1.5, zorder=20)
    
    # Add "Cut Line" text label
    # Place it near the middle of the cut line
    mid_z = height / 2
    cut_x, cut_y, _ = wrap_point(0, mid_z, width, radius, offset)
    ax1.text(cut_x, cut_y - 8, mid_z, "Cut", color='red', fontsize=10, ha='right')
        
    # Set limits - tighter to reduce left whitespace
    margin = 0.5
    ax1.set_xlim(-radius - margin, radius + margin)
    ax1.set_ylim(-radius - margin, radius + margin)
    ax1.set_zlim(-margin, height + margin)
    
    # Adjust view
    ax1.view_init(elev=20, azim=-90)
    
    # Reduce margins around the 3D plot
    ax1.set_box_aspect([1, 1, 1.2])  # Make it more compact
    
    ax1.set_title("Cylindrical Crypt", fontsize=12, y=0.9)

    # --- Right: 2D Flat ---
    ax2 = fig.add_subplot(122)
    ax2.axis('off')
    
    # Draw 2D Polygons
    for p in polygons:
        poly = Polygon(p, facecolor='white', edgecolor='#555555', linewidth=0.5)
        ax2.add_patch(poly)
        
    # Draw Red Cut Lines
    # Left side (Original cut segments)
    for p_start, p_end in cut_segments:
        ax2.plot([p_start[0], p_end[0]], [p_start[1], p_end[1]], color='red', linewidth=1.5)
        
    # Right side (Duplicated at x + width)
    # Since visual wrapping might not be exact per vertex, let's just shift the cut segments by 'width'
    for p_start, p_end in cut_segments:
        ax2.plot([p_start[0] + width, p_end[0] + width], [p_start[1], p_end[1]], color='red', linewidth=1.5, linestyle='--')

    ax2.text(0, -size, "Cut Edge", color='red', ha='center', fontsize=10)
    ax2.text(width, -size, "Cut Edge", color='red', ha='center', fontsize=10)

    ax2.set_aspect('equal')
    ax2.autoscale_view()
    # Add padding
    ax2.set_xlim(-2 * size, width + 2 * size)
    ax2.set_ylim(-size, height + size)
    
    ax2.set_title("Unrolled 2D Lattice", fontsize=12,y=0.92)

    # --- Arrow ---
    # Draw arrow from Right of Left plot to Left of Right plot
    # Coordinates in Figure fraction (0,0 to 1,1)
    
    trans_arrow = FancyArrowPatch(
        (0.42, 0.48), (0.51, 0.48), # Start (x,y), End (x,y)
        transform=fig.transFigure,
        arrowstyle="simple,head_width=10,head_length=10",
        color="black",
        linewidth=1,
        zorder=30
    )
    fig.patches.append(trans_arrow)
    fig.text(0.46, 0.5, "Transform", ha='center', fontsize=10)
    
    plt.tight_layout(pad=0.5)
    plt.savefig('cylinder_transform.png', dpi=300, facecolor='white', bbox_inches='tight', pad_inches=0.05)

if __name__ == "__main__":
    draw_figure()
