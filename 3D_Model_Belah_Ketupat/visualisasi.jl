using PlotlyJS
using JSON

# Function to read the binary model file.
function read_model(filename)
    open(filename, "r") do io
        # Read total vertices (Int32)
        totalVertices = read(io, Int32)
        # Allocate an array for vertex data: totalVertices * 3 Float32 values
        nvertex = totalVertices * 3
        vertexData = Vector{Float32}(undef, nvertex)
        read!(io, vertexData)
        # Reshape into a (totalVertices x 3) matrix (each row is [x y z])
        vertices = reshape(vertexData, (3, totalVertices))'  # transpose so each row is one vertex

        # Read total faces (Int32)
        totalFaces = read(io, Int32)
        # Allocate an array for face data: each face has 4 Int32 values (3 indices and 1 color)
        nface = totalFaces * 4
        faceData = Vector{Int32}(undef, nface)
        read!(io, faceData)
        # Reshape into (totalFaces x 4) matrix
        faces = reshape(faceData, (4, totalFaces))'
        return vertices, faces
    end
end

# Read the model from the binary file.
vertices, faces = read_model("cube_chess_model.bin")

# Separate vertex coordinates.
x = vertices[:, 1]
y = vertices[:, 2]
z = vertices[:, 3]

# Extract face indices and color codes.
# Note: We assume the indices are 0-indexed (as produced in C).
i_idx = faces[:, 1]
j_idx = faces[:, 2]
k_idx = faces[:, 3]
face_color_codes = faces[:, 4]

# Map color codes to color strings.
color_map = Dict(0 => "darkgreen", 1 => "lightgreen", 2 => "orange")
face_colors = [get(color_map, c, "gray") for c in face_color_codes]

# Create a mesh3d trace.
mesh = mesh3d(
    x = x,
    y = y,
    z = z,
    i = i_idx,
    j = j_idx,
    k = k_idx,
    facecolor = face_colors,
    opacity = 1.0,
    flatshading = true
)

plt = Plot(mesh, Layout(title = "Coded by Arief, Idea by Naila", 
    scene = attr(aspectmode = "data")))
display(plt)

# Serialize plot data and layout to JSON strings.
data_str = JSON.json(plt.data)
layout_str = JSON.json(plt.layout)

# Build an HTML document that loads Plotly from the CDN, displays a header title above the plot, and renders the plot.
html_str = """
<html>
<head>
  <meta charset="utf-8" />
  <title>Belah Ketupat</title>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
  <h1 style="text-align:center;">Selamat Hari Raya Idul Fitri, Mohon Maaf Lahir dan Batin</h1>
  <div id="plotly-div" style="width:100%;height:90vh;"></div>
  <script type="text/javascript">
    var plotData = $data_str;
    var plotLayout = $layout_str;
    Plotly.newPlot('plotly-div', plotData, plotLayout);
  </script>
</body>
</html>
"""

# Write the HTML string to a file.
open("selamat_hari_raya.html", "w") do file
    write(file, html_str)
end
