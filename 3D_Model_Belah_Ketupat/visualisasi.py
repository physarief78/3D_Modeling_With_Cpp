import struct
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Nama file biner yang dihasilkan oleh program C
filename = "cube_chess_model.bin"

# Baca seluruh konten file biner
with open(filename, "rb") as f:
    data = f.read()

offset = 0

# Baca jumlah vertex (integer)
total_vertices = struct.unpack("i", data[offset:offset+4])[0]
offset += 4

# Baca array vertex (setiap vertex: 3 float: x, y, z)
vertices = []
for _ in range(total_vertices):
    x, y, z = struct.unpack("fff", data[offset:offset+12])
    vertices.append((x, y, z))
    offset += 12
vertices = np.array(vertices)

# Baca jumlah face (integer)
total_faces = struct.unpack("i", data[offset:offset+4])[0]
offset += 4

# Baca array face (setiap face: 3 integer indeks + 1 integer warna)
faces = []
face_colors = []
for _ in range(total_faces):
    v1, v2, v3, color = struct.unpack("iiii", data[offset:offset+16])
    faces.append((v1, v2, v3))
    face_colors.append(color)
    offset += 16

# Buat list segitiga (list of triangles) untuk plotting
triangles = []
colors = []
for face, col in zip(faces, face_colors):
    triangle = [vertices[face[0]], vertices[face[1]], vertices[face[2]]]
    triangles.append(triangle)
    # Atur warna sesuai atribut: 0 -> darkgreen, 1 -> lightgreen, 2 -> orange (pita)
    if col == 0:
        colors.append("darkgreen")
    elif col == 1:
        colors.append("lightgreen")
    elif col == 2:
        colors.append("orange")
    else:
        colors.append("gray")

# Plot menggunakan matplotlib
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection="3d")

# Buat Poly3DCollection dengan segitiga dan warna masing-masing
mesh = Poly3DCollection(triangles, edgecolor="k", linewidths=0.5, alpha=0.9)
mesh.set_facecolor(colors)
ax.add_collection3d(mesh)

# Atur skala axis berdasarkan data vertex
all_pts = vertices.flatten()
ax.auto_scale_xyz(all_pts, all_pts, all_pts)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_xlim(-1.75, 1.75)
ax.set_ylim(-1.75, 1.75)
plt.title("Nailaa udah siap nih, nanti tinggal di finalisasi :v")

plt.show()
