import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def load_obj(filename):
    """
    Load vertices and faces from an OBJ file.
    Vertex lines: "v x y z"
    Face lines: "f i j k" (1-indexed)
    """
    vertices = []
    faces = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if parts[0] == 'v':
                vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elif parts[0] == 'f':
                faces.append([int(parts[1]) - 1, int(parts[2]) - 1, int(parts[3]) - 1])
    return np.array(vertices), np.array(faces)

def load_magnetic_field_3d(filename):
    """
    Load 3D magnetic field data from a file.
    Expected file format:
      x y z Bx By Bz
    with a header line.
    Returns:
      points: N x 3 array of observation points.
      fields: N x 3 array of corresponding magnetic field vectors.
    """
    data = np.loadtxt(filename, skiprows=1)
    points = data[:, :3]
    fields = data[:, 3:]
    return points, fields

def plot_mesh(ax, vertices, faces, color):
    """
    Plot a mesh on the provided Axes3D object using a Poly3DCollection.
    """
    polys = []
    for face in faces:
        poly = vertices[face]
        polys.append(poly)
    mesh_collection = Poly3DCollection(polys, facecolor=color, edgecolor='k', alpha=0.7)
    ax.add_collection3d(mesh_collection)

def main():
    # Load the torus OBJ files.
    vertices1, faces1 = load_obj("torus1.obj")
    vertices2, faces2 = load_obj("torus2.obj")
    
    # Load the 3D magnetic field data.
    points, fields = load_magnetic_field_3d("magnetic_field_3d.txt")
    
    # Create a figure with two subplots: left is 3D, right is a 2D surface plot.
    fig = plt.figure(figsize=(16, 8))
    ax3d = fig.add_subplot(121, projection='3d')
    ax2d = fig.add_subplot(122)
    
    # --- 3D Visualization ---
    # Plot the torus pair.
    plot_mesh(ax3d, vertices1, faces1, color='red')
    plot_mesh(ax3d, vertices2, faces2, color='blue')
    
    # Sample and plot magnetic field vectors in 3D.
    sample_rate = 2  # adjust to plot fewer or more arrows
    sample_indices = np.arange(0, points.shape[0], sample_rate)
    sample_points = points[sample_indices]
    sample_fields = fields[sample_indices]
    ax3d.quiver(sample_points[:, 0], sample_points[:, 1], sample_points[:, 2],
                 sample_fields[:, 0], sample_fields[:, 1], sample_fields[:, 2],
                 length=0.25, normalize=True, color='green')
    
    # Set axis labels and limits for 3D plot.
    all_data = np.concatenate([vertices1, vertices2, points], axis=0)
    ax3d.set_xlim(all_data[:,0].min(), all_data[:,0].max())
    ax3d.set_ylim(all_data[:,1].min(), all_data[:,1].max())
    ax3d.set_zlim(all_data[:,2].min(), all_data[:,2].max())
    ax3d.set_xlabel("X")
    ax3d.set_ylabel("Y")
    ax3d.set_zlabel("Z")
    ax3d.set_title("3D Helmholtz's Coil Simulation")
    
    # --- 2D Surface Magnetic Field Plot ---
    # Extract a horizontal slice at z ≈ 0.
    unique_z = np.unique(points[:, 2])
    z_target = unique_z[np.argmin(np.abs(unique_z))]
    tol = 1e-3
    slice_mask = np.abs(points[:, 2] - z_target) < tol
    slice_points = points[slice_mask]
    slice_fields = fields[slice_mask]
    
    # Compute the magnetic field magnitude at each observation point.
    field_magnitude = np.sqrt(np.sum(slice_fields**2, axis=1))
    
    # For a structured grid, sort the slice data by x then y.
    sort_idx = np.lexsort((slice_points[:, 1], slice_points[:, 0]))
    slice_points = slice_points[sort_idx]
    field_magnitude = field_magnitude[sort_idx]
    slice_fields = slice_fields[sort_idx]
    
    # Extract unique x and y values.
    unique_x = np.unique(slice_points[:, 0])
    unique_y = np.unique(slice_points[:, 1])
    nx = unique_x.size
    ny = unique_y.size
    
    # Reshape the field magnitude into a 2D array.
    Z_field = field_magnitude.reshape((nx, ny)).T
    
    # Reshape the horizontal (x and y) components of the magnetic field.
    U_field = slice_fields[:, 0].reshape((nx, ny)).T
    V_field = slice_fields[:, 1].reshape((nx, ny)).T
    
    # Create a meshgrid for the x-y plane.
    X, Y = np.meshgrid(unique_x, unique_y)
    
    # Plot the 2D contour/fill plot.
    contour = ax2d.contourf(X, Y, Z_field, 50, cmap='RdBu')
    fig.colorbar(contour, ax=ax2d, label="Magnetic Field Magnitude (T)")
    ax2d.set_xlabel("X")
    ax2d.set_ylabel("Y")
    ax2d.set_title("2D Magnetic Field at z = {:.3f}".format(z_target))
    
    # Overlay vector field on the 2D surface.
    # Adjust the sampling (skip) for clarity.
    skip = slice(None, None, 2)  # change step value as needed
    ax2d.quiver(X[skip, skip], Y[skip, skip],
                U_field[skip, skip], V_field[skip, skip],
                color='white', scale=50, width=0.003)
    
    # Save the figure to a file before showing.
    plt.savefig("helmholtz_field.png", dpi=300, bbox_inches='tight')
    print("Figure saved as helmholtz_field.png")
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
