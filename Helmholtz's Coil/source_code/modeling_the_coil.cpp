#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// -- Geometry generation for torus meshes --

struct Vertex {
    float x, y, z;
};

struct Mesh {
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
};

Mesh createTorus(float R, float r, int numSides, int numRings) {
    Mesh mesh;

    for (int i = 0; i < numRings; ++i) {
        float u = i * 2.0f * M_PI / numRings;
        float cosU = std::cos(u);
        float sinU = std::sin(u);

        for (int j = 0; j < numSides; ++j) {
            float v = j * 2.0f * M_PI / numSides;
            float cosV = std::cos(v);
            float sinV = std::sin(v);

            Vertex vert;
            // Parameterization: a torus with tube along y.
            vert.x = (R + r * cosV) * cosU;
            vert.y = r * sinV;
            vert.z = (R + r * cosV) * sinU;
            mesh.vertices.push_back(vert);
        }
    }

    // Generate triangle faces from the grid
    for (int i = 0; i < numRings; ++i) {
        int nextRing = (i + 1) % numRings;
        for (int j = 0; j < numSides; ++j) {
            int nextSide = (j + 1) % numSides;
            int current    = i * numSides + j;
            int right      = i * numSides + nextSide;
            int below      = nextRing * numSides + j;
            int belowRight = nextRing * numSides + nextSide;

            mesh.indices.push_back(current);
            mesh.indices.push_back(below);
            mesh.indices.push_back(right);

            mesh.indices.push_back(right);
            mesh.indices.push_back(below);
            mesh.indices.push_back(belowRight);
        }
    }
    return mesh;
}

// Rotate a mesh about the z-axis by a given angle (in radians)
void rotateMeshZ(Mesh &mesh, float angle) {
    float cosA = std::cos(angle);
    float sinA = std::sin(angle);
    for (auto &vertex : mesh.vertices) {
        float x = vertex.x;
        float y = vertex.y;
        vertex.x = x * cosA - y * sinA;
        vertex.y = x * sinA + y * cosA;
        // vertex.z remains unchanged.
    }
}

// Save the mesh to an OBJ file.
void saveMeshToOBJ(const Mesh& mesh, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
        return;
    }
    for (const auto& v : mesh.vertices) {
        file << "v " << v.x << " " << v.y << " " << v.z << "\n";
    }
    for (size_t i = 0; i < mesh.indices.size(); i += 3) {
        file << "f " << mesh.indices[i] + 1 << " " 
             << mesh.indices[i + 1] + 1 << " " 
             << mesh.indices[i + 2] + 1 << "\n";
    }
    file.close();
    std::cout << "Saved mesh to " << filename << std::endl;
}

// -- Magnetic field simulation for Helmholtz coils --

// A simple 3D vector for double-precision arithmetic.
struct Vec3 {
    double x, y, z;
};

Vec3 add(const Vec3 &a, const Vec3 &b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3 sub(const Vec3 &a, const Vec3 &b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 cross(const Vec3 &a, const Vec3 &b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}
double norm(const Vec3 &a) {
    return std::sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

// Compute the magnetic field at a point due to a circular current loop.
// The coil is approximated by 'segments' straight current elements.
Vec3 computeCoilField(const Vec3 &center, double radius, double I, int segments, const Vec3 &point) {
    const double mu0 = 4 * M_PI * 1e-7;  // permeability of free space
    double dtheta = 2 * M_PI / segments;
    Vec3 B = {0, 0, 0};
    for (int i = 0; i < segments; ++i) {
        double theta = i * dtheta;
        double theta_next = (i + 1) * dtheta;
        // The coil is assumed to lie in the y-z plane (with constant x = center.x)
        Vec3 pos = { center.x, center.y + radius * std::cos(theta), center.z + radius * std::sin(theta) };
        Vec3 pos_next = { center.x, center.y + radius * std::cos(theta_next), center.z + radius * std::sin(theta_next) };
        // Current element dl
        Vec3 dl = { pos_next.x - pos.x, pos_next.y - pos.y, pos_next.z - pos.z };
        // Vector from the current element to the observation point.
        Vec3 r_vec = { point.x - pos.x, point.y - pos.y, point.z - pos.z };
        double r = norm(r_vec);
        if (r < 1e-6) continue;  // avoid singularity
        Vec3 dB = cross(dl, r_vec);
        double factor = (mu0 * I) / (4 * M_PI * std::pow(r, 3));
        dB.x *= factor;
        dB.y *= factor;
        dB.z *= factor;
        B = add(B, dB);
    }
    return B;
}

int main() {
    // -- Geometry parameters --
    float majorRadius = 2.0f;
    float minorRadius = 0.5f;
    int numSides = 30;
    int numRings = 30;

    // Create the torus meshes
    Mesh torus1 = createTorus(majorRadius, minorRadius, numSides, numRings);
    Mesh torus2 = createTorus(majorRadius, minorRadius, numSides, numRings);

    // Rotate torus1 by -90° about z so that its “hole” points toward positive x.
    rotateMeshZ(torus1, -M_PI / 2);
    // Rotate torus2 by +90° about z so that its “hole” points toward negative x.
    rotateMeshZ(torus2, M_PI / 2);

    // For Helmholtz configuration, the two coils (approximated by the torus centerlines)
    // are translated along the x-axis. Ideally, the separation between coil centers equals
    // the effective coil radius. Here we use:
    double effectiveCoilRadius = majorRadius + minorRadius; // effective radius of the current loop
    double translationDistance = 2;    // Helmholtz: separation = coil radius
    for (auto &vertex : torus1.vertices) {
        vertex.x -= translationDistance;  // Move torus1 to the left.
    }
    for (auto &vertex : torus2.vertices) {
        vertex.x += translationDistance;  // Move torus2 to the right.
    }

    // Save the torus meshes (optional for visualization)
    saveMeshToOBJ(torus1, "torus1.obj");
    saveMeshToOBJ(torus2, "torus2.obj");

    // -- Helmholtz Coil Magnetic Field Simulation (1D) --
    // Approximate each torus as a circular current loop.
    // Coil1 is centered at (-translationDistance, 0, 0), Coil2 at (translationDistance, 0, 0).
    Vec3 coil1_center = { -2, 0, 0 };
    Vec3 coil2_center = {  2, 0, 0 };

    // Set the current in each coil (in Amperes)
    double I = 5.0;
    int segments = 1000;  // resolution for numerical integration

    // (Optional) Compute and save magnetic field along the x-axis.
    std::vector<Vec3> observationPoints1D;
    int numPoints1D = 30;
    for (int i = 0; i < numPoints1D; ++i) {
        double x = -effectiveCoilRadius + i * (2 * effectiveCoilRadius) / (numPoints1D - 1);
        observationPoints1D.push_back({x, 0, 0});
    }
    std::ofstream outFile1D("magnetic_field.txt");
    if (!outFile1D) {
        std::cerr << "Error: Could not open magnetic_field.txt for writing." << std::endl;
        return 1;
    }
    outFile1D << "Magnetic field at observation points along x-axis (in Tesla):\n";
    for (const auto &pt : observationPoints1D) {
        Vec3 B1 = computeCoilField(coil1_center, effectiveCoilRadius, I, segments, pt);
        Vec3 B2 = computeCoilField(coil2_center, effectiveCoilRadius, I, segments, pt);
        Vec3 B_total = { B1.x + B2.x, B1.y + B2.y, B1.z + B2.z };
        outFile1D << "At (" << pt.x << ", " << pt.y << ", " << pt.z << "): "
                  << "B = (" << B_total.x << ", " << B_total.y << ", " << B_total.z << ") T\n";
    }
    outFile1D.close();
    std::cout << "1D magnetic field data saved to magnetic_field.txt" << std::endl;

    // -- Helmholtz Coil Magnetic Field Simulation (3D) --
    // Create a 3D grid of observation points.
    double grid_min = -effectiveCoilRadius * 3.0;
    double grid_max =  effectiveCoilRadius * 3.0;
    int grid_points = 22; // number of points along each axis
    double step = (grid_max - grid_min) / (grid_points - 1);

    std::ofstream outFile3D("magnetic_field_3d.txt");
    if (!outFile3D) {
        std::cerr << "Error: Could not open magnetic_field_3d.txt for writing." << std::endl;
        return 1;
    }
    outFile3D << "x y z Bx By Bz\n";
    for (int i = 0; i < grid_points; ++i) {
        double x = grid_min + i * step;
        for (int j = 0; j < grid_points; ++j) {
            double y = grid_min + j * step;
            for (int k = 0; k < grid_points; ++k) {
                double z = grid_min + k * step;
                Vec3 pt = {x, y, z};
                Vec3 B1 = computeCoilField(coil1_center, effectiveCoilRadius, I, segments, pt);
                Vec3 B2 = computeCoilField(coil2_center, effectiveCoilRadius, I, segments, pt);
                Vec3 B_total = { B1.x + B2.x, B1.y + B2.y, B1.z + B2.z };
                outFile3D << x << " " << y << " " << z << " "
                          << B_total.x << " " << B_total.y << " " << B_total.z << "\n";
            }
        }
    }
    outFile3D.close();
    std::cout << "3D magnetic field data saved to magnetic_field_3d.txt" << std::endl;

    return 0;
}
