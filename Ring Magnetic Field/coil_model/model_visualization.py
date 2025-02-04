import numpy as np
import matplotlib.pyplot as plt

# Function to load coil geometry from the file
def load_coil_geometry(filename):
    coil1 = []
    coil2 = []
    # coil3 = []
    # coil4 = []
    
    with open(filename, 'r') as file:
        lines = file.readlines()
        current_coil = None
        
        for line in lines:
            line = line.strip()
            if line.startswith("Coil 1 Points:"):
                current_coil = coil1
            elif line.startswith("Coil 2 Points:"):
                current_coil = coil2
            # elif line.startswith("Coil 3 Points:"):
            #     current_coil = coil3
            # elif line.startswith("Coil 4 Points:"):
            #     current_coil = coil4
            elif line and current_coil is not None:
                x, y, z = map(float, line.split(", "))
                current_coil.append([x, y, z])
    
    return np.array(coil1), np.array(coil2)

# Function to visualize the coil geometry
def visualize_coils(coil1, coil2):

    # 3D plot for coil geometry
    fig = plt.figure(figsize=(14, 6))

    ax = fig.add_subplot(111, projection='3d')
    ax.plot(coil1[:, 0], coil1[:, 1], coil1[:, 2], label="Coil 1", color='red')
    ax.plot(coil2[:, 0], coil2[:, 1], coil2[:, 2], label="Coil 2", color='blue')
    # ax.plot(coil3[:, 0], coil3[:, 1], coil3[:, 2], label="Coil 3", color='green')
    # ax.plot(coil4[:, 0], coil4[:, 1], coil4[:, 2], label="Coil 4", color='purple')

    ax.set_title("3D Coil Geometry")
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.set_zlim(-10, 10)
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")
    ax.legend()

    plt.tight_layout()
    plt.show()

# Main function
def main():
    coil_geometry_file = "serial_maxwell_coil_model.txt"

    # Load coil geometry data
    coil1, coil2 = load_coil_geometry(coil_geometry_file)

    # Visualize the coils
    visualize_coils(coil1, coil2)

if __name__ == "__main__":
    main()
