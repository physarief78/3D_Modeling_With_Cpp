import numpy as np
import matplotlib.pyplot as plt

def load_coil_geometry(filename):
    coils = {1: [], 2: [], 3: [], 4: []}
    current_coil = None
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("Coil"):
                current_coil = int(line.split()[1])
            elif line:
                coils[current_coil].append(list(map(float, line.split(", "))))
    return [np.array(coils[i]) for i in range(1, 5)]

def load_field_data(filename):
    return np.array([list(map(float, line.strip().split(','))) for line in open(filename)])

def visualize_coils_and_field(coils, field_data, save_image=False):
    field_points, field_vectors = field_data[:, :3], field_data[:, 3:]
    magnitude = np.linalg.norm(field_vectors, axis=1)
    x, y = field_points[:, 0], field_points[:, 1]
    X, Y = np.meshgrid(np.unique(x), np.unique(y))
    Z = np.zeros_like(X)
    for i, (xi, yi) in enumerate(zip(x, y)):
        Z[np.where(np.unique(y) == yi)[0][0], np.where(np.unique(x) == xi)[0][0]] = magnitude[i]

    fig = plt.figure(figsize=(14, 6))
    ax1 = fig.add_subplot(121, projection='3d')
    colors = ['red', 'blue', 'green', 'orange']
    for i, coil in enumerate(coils):
        ax1.plot(coil[:, 0], coil[:, 1], coil[:, 2], label=f"Coil {i+1}", color=colors[i])
    ax1.quiver(field_points[:, 0], field_points[:, 1], field_points[:, 2], field_vectors[:, 0], field_vectors[:, 1], field_vectors[:, 2],
               length=0.75, normalize=True, color='magenta', label="Magnetic Field")
    ax1.set_title("3D Coil Geometry and Magnetic Field Vector")
    ax1.set_xlabel("X-axis")
    ax1.set_ylabel("Y-axis")
    ax1.set_zlabel("Z-axis")
    ax1.legend()

    ax2 = fig.add_subplot(122)
    surface = ax2.contourf(X, Y, Z, levels=50, cmap='inferno')
    plt.colorbar(surface, ax=ax2, label="Magnetic Field Magnitude")
    ax2.set_title("2D Magnetic Field Magnitude Surface")
    ax2.set_xlabel("X-axis")
    ax2.set_ylabel("Y-axis")

    plt.tight_layout()
    if save_image:
        plt.savefig("two_coil_pair_and_field_visualization.png", dpi=300, bbox_inches='tight')
    plt.show()

def main():
    coil1, coil2, coil3, coil4 = load_coil_geometry("two_coil_pair_geometry.txt")
    field_data = load_field_data("two_coil_pair_magnetic_field_data.txt")
    visualize_coils_and_field([coil1, coil2, coil3, coil4], field_data, save_image=True)

if __name__ == "__main__":
    main()