#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

// Structure to hold a vertex. We also store the spherical angles.
struct Vertex {
    float x, y, z;
    float theta, phi; // Spherical coordinates (polar and azimuthal angles)
};

// Structure to hold the mesh data.
struct Mesh {
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
};

// Function to create the base sphere mesh (storing theta and phi for later use).
Mesh createSphereMesh(float radius, unsigned int stacks, unsigned int slices) {
    Mesh mesh;
    for (unsigned int i = 0; i <= stacks; i++) {
        float theta = i * M_PI / stacks;
        float sinTheta = std::sin(theta);
        float cosTheta = std::cos(theta);
        for (unsigned int j = 0; j <= slices; j++) {
            float phi = j * 2 * M_PI / slices;
            Vertex vertex;
            vertex.theta = theta;
            vertex.phi = phi;
            // Compute the initial Cartesian coordinates.
            vertex.x = radius * sinTheta * std::cos(phi);
            vertex.y = radius * sinTheta * std::sin(phi);
            vertex.z = radius * cosTheta;
            mesh.vertices.push_back(vertex);
        }
    }
    // Generate indices for the triangles (connectivity remains constant).
    for (unsigned int i = 0; i < stacks; i++) {
        for (unsigned int j = 0; j < slices; j++) {
            unsigned int first = i * (slices + 1) + j;
            unsigned int second = first + slices + 1;
            // First triangle.
            mesh.indices.push_back(first);
            mesh.indices.push_back(second);
            mesh.indices.push_back(first + 1);
            // Second triangle.
            mesh.indices.push_back(second);
            mesh.indices.push_back(second + 1);
            mesh.indices.push_back(first + 1);
        }
    }
    return mesh;
}

// Eigenstate function that returns a radial displacement for given angles and time.
// The deformation is defined as:
//   Δr(θ,φ,t) = A * sin(eigenMode * θ) * cos((eigenMode - 1)*φ + t)
// By adjusting 'eigenMode', you change the spatial frequency (or "level") of the eigenstate.
float eigenstate(float theta, float phi, float time, float amplitude, int eigenMode) {
    return amplitude * std::sin(eigenMode * theta) * std::cos((eigenMode - 1) * phi + time);
}

// Update the Cartesian positions of the mesh vertices based on the eigenstate.
void updateVertexPositions(Mesh &mesh, float baseRadius, float amplitude, float time, int eigenMode) {
    for (auto &v : mesh.vertices) {
        float displacement = eigenstate(v.theta, v.phi, time, amplitude, eigenMode);
        float newRadius = baseRadius + displacement;
        v.x = newRadius * std::sin(v.theta) * std::cos(v.phi);
        v.y = newRadius * std::sin(v.theta) * std::sin(v.phi);
        v.z = newRadius * std::cos(v.theta);
    }
}

int main() {
    // Parameters.
    float baseRadius = 1.0f;
    float amplitude = 0.2f;         // Amplitude of the eigenstate deformation.
    unsigned int stacks = 30;       // Latitudinal subdivisions.
    unsigned int slices = 30;       // Longitudinal subdivisions.
    float timeStart = 0.0f;
    float timeEnd = 8 * M_PI;       // One full cycle.
    unsigned int numFrames = 500;   // Total number of frames.
    float dt = (timeEnd - timeStart) / numFrames;

    // Set the eigenMode level; e.g., eigenMode = 3 gives the original deformation pattern.
    int eigenMode = 5;  // You can adjust this value to change the oscillation frequency.

    // Create the base sphere mesh.
    Mesh sphereMesh = createSphereMesh(baseRadius, stacks, slices);

    // Open a file to write JSON data.
    std::ofstream outFile("animation_data.json");
    if (!outFile) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return 1;
    }

    // Write JSON header with fixed floating-point format.
    outFile << std::fixed << std::setprecision(6);
    outFile << "{\n";

    // Write the constant indices (mesh connectivity).
    outFile << "  \"indices\": [\n    ";
    for (size_t i = 0; i < sphereMesh.indices.size(); i++) {
        outFile << sphereMesh.indices[i];
        if (i < sphereMesh.indices.size() - 1)
            outFile << ", ";
        if ((i+1) % 10 == 0)
            outFile << "\n    ";
    }
    outFile << "\n  ],\n";

    // Write frame data.
    outFile << "  \"frames\": [\n";
    for (unsigned int frame = 0; frame <= numFrames; frame++) {
        float currentTime = timeStart + frame * dt;
        // Update vertex positions for this frame.
        updateVertexPositions(sphereMesh, baseRadius, amplitude, currentTime, eigenMode);

        outFile << "    {\n";
        outFile << "      \"time\": " << currentTime << ",\n";
        outFile << "      \"vertices\": [\n        ";
        for (size_t i = 0; i < sphereMesh.vertices.size(); i++) {
            const auto &v = sphereMesh.vertices[i];
            outFile << "[" << v.x << ", " << v.y << ", " << v.z << "]";
            if (i < sphereMesh.vertices.size() - 1)
                outFile << ", ";
            if ((i+1) % 5 == 0)
                outFile << "\n        ";
        }
        outFile << "\n      ]\n    }";
        if (frame < numFrames)
            outFile << ",";
        outFile << "\n";
    }
    outFile << "  ]\n";

    outFile << "}\n";
    outFile.close();

    std::cout << "Animation data saved to animation_data.json" << std::endl;
    return 0;
}
