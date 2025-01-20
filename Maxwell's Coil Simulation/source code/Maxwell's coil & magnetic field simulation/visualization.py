import numpy as np
import matplotlib.pyplot as plt

# Function to load coil geometry from the file
def load_coil_geometry(filename):
    coil1 = []
    coil2 = []
    coil3 = []
    coil4 = []
    
    with open(filename, 'r') as file:
        lines = file.readlines()
        current_coil = None
        
        for line in lines:
            line = line.strip()
            if line.startswith("Coil 1 Points:"):
                current_coil = coil1
            elif line.startswith("Coil 2 Points:"):
                current_coil = coil2
            elif line.startswith("Coil 3 Points:"):
                current_coil = coil3
            elif line.startswith("Coil 4 Points:"):
                current_coil = coil4
            elif line and current_coil is not None:
                x, y, z = map(float, line.split(", "))
                current_coil.append([x, y, z])
    
    return np.array(coil1), np.array(coil2), np.array(coil3), np.array(coil4)

# Function to load magnetic field data from the file
def load_field_data(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            values = list(map(float, line.strip().split(',')))
            data.append(values)
    return np.array(data)

# Function to visualize the coil geometry
def visualize_coils(coil1, coil2, coil3, coil4, field_data):
    #Extract field components
    field_points = field_data[:, :3]
    field_vectors = field_data[:, 3:]

    # 3D plot for coil geometry
    fig = plt.figure(figsize=(14, 6))

    ax1 = fig.add_subplot(121, projection='3d')
    ax1.plot(coil1[:, 0], coil1[:, 1], coil1[:, 2], label="Coil 1", color='red')
    ax1.plot(coil2[:, 0], coil2[:, 1], coil2[:, 2], label="Coil 2", color='blue')
    ax1.plot(coil3[:, 0], coil3[:, 1], coil3[:, 2], label="Coil 3", color='green')
    ax1.plot(coil4[:, 0], coil4[:, 1], coil4[:, 2], label="Coil 4", color='purple')

    ax1.quiver(
        field_points[:, 0], field_points[:, 1], field_points[:, 2],
        field_vectors[:, 0], field_vectors[:, 1], field_vectors[:, 2],
        length=0.3, normalize=True, color='magenta', label="Magnetic Field"
    )

    ax1.set_title("3D Coil Geometry")
    ax1.set_xlim(-7, 7)
    ax1.set_ylim(-7, 7)
    ax1.set_zlim(-7, 7)
    ax1.set_xlabel("X-axis")
    ax1.set_ylabel("Y-axis")
    ax1.set_zlabel("Z-axis")
    ax1.legend()

    # 2D surface plot for magnetic field magnitude
    magnitude = np.linalg.norm(field_vectors, axis=1)
    x = field_points[:, 0]
    #y = field_points[:, 1]
    z = field_points[:, 2]

    # Reshape for 2D grid plotting
    unique_x = np.unique(x)
    #unique_y = np.unique(y)
    unique_z = np.unique(z)
    #X, Y = np.meshgrid(unique_x, unique_y)
    X, Z = np.meshgrid(unique_x, unique_z)
    #Z = np.zeros_like(X)
    Y = np.zeros_like(X)

    for i in range(len(x)):
        xi = np.where(unique_x == x[i])[0][0]
        #yi = np.where(unique_y == y[i])[0][0]
        zi = np.where(unique_z == z[i])[0][0]
        #Z[yi, xi] = magnitude[i]
        Y[zi, xi] = magnitude[i]

    ax2 = fig.add_subplot(122)
    #surface = ax2.contourf(X, Y, Z, levels=50, cmap='inferno')
    surface = ax2.contourf(X, Z, Y, levels=100, cmap='inferno')
    plt.colorbar(surface, ax=ax2, label="Magnetic Field Magnitude")

    ax2.set_title("2D Magnetic Field Magnitude Surface")
    ax2.set_xlim(field_points[:, 0].min(), field_points[:, 0].max())
    ax2.set_ylim(field_points[:, 2].min(), field_points[:, 2].max())
    ax2.set_xlabel("X-axis")
    ax2.set_ylabel("Z-axis")

    plt.tight_layout()
    plt.show()

# Main function
def main():
    coil_geometry_file = "4_maxwell_coil_model.txt"
    field_data_file = "4_coil_magnetic_field_data.txt"

    # Load coil geometry data
    coil1, coil2, coil3, coil4 = load_coil_geometry(coil_geometry_file)
    field_data = load_field_data(field_data_file)

    # Visualize the coils
    visualize_coils(coil1, coil2, coil3, coil4, field_data)

if __name__ == "__main__":
    main()
