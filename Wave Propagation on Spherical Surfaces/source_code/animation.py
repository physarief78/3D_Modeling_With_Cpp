import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.animation as animation

# Load the animation data from the JSON file.
with open("animation_data.json", "r") as f:
    data = json.load(f)

indices = data["indices"]
frames = data["frames"]

# Convert the flat indices list into a list of triangles.
# Each triangle is a triplet of indices.
triangles_indices = [indices[i:i+3] for i in range(0, len(indices), 3)]

def get_triangles_for_frame(frame):
    """Constructs a list of triangles (each triangle is a list of 3 vertices)
    for the given frame data."""
    vertices = frame["vertices"]
    triangles = []
    for tri in triangles_indices:
        # Each vertex is a list [x, y, z]
        triangle = [vertices[idx] for idx in tri]
        triangles.append(triangle)
    return triangles

# Set up the figure and 3D axes.
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set axis limits (adjust based on your sphere size and deformation amplitude).
ax.set_xlim(-7.5, 7.5)
ax.set_ylim(-7.5, 7.5)
ax.set_zlim(-7.5, 7.5)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

# Adjustable view parameters.
initial_elev = 20    # Initial elevation angle in degrees.
initial_azim = 45   # Initial azimuth angle in degrees.
rotation_speed = 0.5   # Degrees to rotate per frame.

# Set the initial view.
ax.view_init(elev=initial_elev, azim=initial_azim)

# Create an initial Poly3DCollection for the first frame.
triangles0 = get_triangles_for_frame(frames[0])
mesh_collection = Poly3DCollection(triangles0, facecolor='cyan', edgecolor='grey', alpha=0.8)
ax.add_collection3d(mesh_collection)

def update(frame_index):
    """Animation update function.
    For each frame, update the mesh vertices and camera view angle."""
    frame = frames[frame_index]
    triangles = get_triangles_for_frame(frame)
    mesh_collection.set_verts(triangles)
    
    # Update the title with current time.
    ax.set_title(f"Wave Propagation in Spherical Surfaces at Time: {frame['time']:.2f}")
    
    # Rotate the camera view by updating the azimuth angle.
    #new_azim = initial_azim + rotation_speed * frame_index
    #ax.view_init(elev=initial_elev, azim=new_azim)
    
    return mesh_collection,

# Create the animation: interval in milliseconds (adjust as needed).
ani = animation.FuncAnimation(fig, update, frames=len(frames), interval=50, blit=False)

# Save the animation to a file.
# Ensure that ffmpeg is installed on your system.
ani.save("wave_propagation_animation.gif", writer="pillow", fps=30)

plt.show()
