import matplotlib.pyplot as plt
import numpy as np

def hexagon(center, size):
    """Generate the vertices of a hexagon given its center and size."""
    angles = np.linspace(0, 2 * np.pi, 7)
    x_hex = center[0] + size * np.cos(angles)
    y_hex = center[1] + size * np.sin(angles)
    return x_hex, y_hex

def lozenge(center, size):
    """Generate the vertices of the bottom right lozenge of the central hexagon."""
    # Bottom right lozenge vertices
    vertices_x = [
        center[0], 
        center[0] + size * np.cos(np.pi / 6),
        center[0] + size * np.cos(np.pi / 6) - size * np.sin(np.pi / 6),
        center[0] - size * np.sin(np.pi / 6)
    ]
    vertices_y = [
        center[1], 
        center[1] + size * np.sin(np.pi / 6),
        center[1] + size * np.sin(np.pi / 6) + size * np.cos(np.pi / 6),
        center[1] + size * np.cos(np.pi / 6)
    ]
    return vertices_x, vertices_y

def plot_motif():
    # Define the size of the hexagon
    size = 1

    # Define the centers of the hexagons
    centers = [
        (0, 0),  # center
        (np.sqrt(3) * size, 0),  # right
        (np.sqrt(3) / 2 * size, -1.5 * size),  # bottom right
    ]

    fig, ax = plt.subplots()

    # Plot the bottom right lozenge of the central hexagon
    x_lozenge, y_lozenge = lozenge(centers[0], size)
    ax.fill(x_lozenge, y_lozenge, edgecolor='black', fill=True, color='lightblue')

    # Plot the two surrounding hexagons
    for center in centers[1:]:
        x_hex, y_hex = hexagon(center, size)
        ax.fill(x_hex, y_hex, edgecolor='black', fill=False)

    # Plot the central hexagon outline
    x_hex_center, y_hex_center = hexagon(centers[0], size)
    ax.plot(x_hex_center, y_hex_center, 'black')

    # Set equal scaling
    ax.set_aspect('equal')

    # Remove axes
    ax.axis('off')

    # Show the plot
    plt.show()
    plt.savefig('hexagon.png')
    plt.close()

plot_motif()
