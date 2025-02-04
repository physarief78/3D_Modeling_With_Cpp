#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// Structure to hold 3D point data
struct Point3D {
    double x, y, z;
};

// Function to generate points for a torus-shaped coil
std::vector<Point3D> generateTorus(double majorRadius, double minorRadius, double centerX, double centerY, double centerZ, int numMajorPoints, int numMinorPoints) {
    std::vector<Point3D> coilPoints;
    double majorAngleStep = 2 * M_PI / numMajorPoints;
    double minorAngleStep = 2 * M_PI / numMinorPoints;

    for (int i = 0; i < numMajorPoints; ++i) {
        double majorAngle = i * majorAngleStep;

        // Center of the minor circle in the torus
        double baseX = centerX + majorRadius * cos(majorAngle);
        double baseY = centerY;
        double baseZ = centerZ + majorRadius * sin(majorAngle);

        for (int j = 0; j < numMinorPoints; ++j) {
            double minorAngle = j * minorAngleStep;

            // Points distributed around the minor circle
            double offsetX = minorRadius * sin(minorAngle) * cos(majorAngle);
            double offsetY = minorRadius * cos(minorAngle);
            double offsetZ = minorRadius * sin(minorAngle) * sin(majorAngle);

            coilPoints.push_back({baseX + offsetX, baseY + offsetY, baseZ + offsetZ});
        }
    }

    return coilPoints;
}

// Function to save coil geometry to a file
void saveCoilGeometryToFile(const std::string& filename, const std::vector<Point3D>& coil1, const std::vector<Point3D>& coil2) {
    std::ofstream outFile(filename);

    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    outFile << "Coil 1 Points:\n";
    for (const auto& point : coil1) {
        outFile << point.x << ", " << point.y << ", " << point.z << "\n";
    }

    outFile << "\nCoil 2 Points:\n";
    for (const auto& point : coil2) {
        outFile << point.x << ", " << point.y << ", " << point.z << "\n";
    }

    outFile.close();
    std::cout << "Coil geometry saved to " << filename << std::endl;
}

// Function to calculate the magnetic field at a point using Biot-Savart's Law
struct Vector3D {
    double x, y, z;

    // Vector addition
    Vector3D operator+(const Vector3D& v) const {
        return {x + v.x, y + v.y, z + v.z};
    }

    // Scalar multiplication
    Vector3D operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }
};

Vector3D calculateMagneticField(const Point3D& p, const std::vector<Point3D>& coil, double current) {
    const double mu0 = 4 * M_PI * 1e-7; // Permeability of free space
    Vector3D B = {0.0, 0.0, 0.0};

    for (size_t i = 0; i < coil.size(); ++i) {
        size_t next = (i + 1) % coil.size();

        // Current segment (dl)
        double dlx = coil[next].x - coil[i].x;
        double dly = coil[next].y - coil[i].y;
        double dlz = coil[next].z - coil[i].z;

        // Vector from the segment to the point p (r)
        double rx = p.x - coil[i].x;
        double ry = p.y - coil[i].y;
        double rz = p.z - coil[i].z;

        double r_mag = std::sqrt(rx * rx + ry * ry + rz * rz);
        double r_mag3 = r_mag * r_mag * r_mag;

        // Biot-Savart law: dB = (mu0 / (4 * pi)) * (I * dl x r) / r^3
        double factor = (mu0 * current) / (4 * M_PI * r_mag3);

        Vector3D dB = {
            factor * (dly * rz - dlz * ry),
            factor * (dlz * rx - dlx * rz),
            factor * (dlx * ry - dly * rx)
        };

        B = B + dB;
    }

    return B;
}

// Function to save magnetic field data to a file
void saveFieldToFile(const std::string& filename, const std::vector<std::pair<Point3D, Vector3D>>& fieldData) {
    std::ofstream outFile(filename);

    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    for (const auto& data : fieldData) {
        outFile << data.first.x << ", " << data.first.y << ", " << data.first.z << ", "
                << data.second.x << ", " << data.second.y << ", " << data.second.z << "\n";
    }

    outFile.close();
    std::cout << "Magnetic field data saved to " << filename << std::endl;
}

int main() {
     // Parameters for Maxwell's coil pair
    double majorRadius = 3.0;    // Major radius of the torus
    double minorRadius = 0.1;    // Minor radius (thickness of the torus)
    double separation = 12.0;     // Distance between the two torus coils
    int numMajorPoints = 200;    // Points along the major circle
    int numMinorPoints = 200;     // Points along the minor circle
    double currentCoil1 = 1.0;   // Current through the first coil (in amperes)
    double currentCoil2 = -1.0;  // Current through the second coil (in amperes, opposite direction)

    // Generate points for the first torus
    std::vector<Point3D> coil1 = generateTorus(majorRadius, minorRadius, separation / 2.0, 0.0, 0.0, numMajorPoints, numMinorPoints);

    // Generate points for the second torus
    std::vector<Point3D> coil2 = generateTorus(majorRadius, minorRadius, -separation / 2.0, 0.0, 0.0, numMajorPoints, numMinorPoints);

    // Save coil geometry
    saveCoilGeometryToFile("a_coil_pair_geometry.txt", coil1, coil2);

    // Define a grid of points to calculate the magnetic field
    std::vector<std::pair<Point3D, Vector3D>> fieldData;
    double gridRange = 10.0; // Range of the grid (from -gridRange to +gridRange)
    double gridStep = 2.0; // Step size for the grid

    for (double x = -gridRange; x <= gridRange; x += gridStep) {
        for (double y = -gridRange; y <= gridRange; y += gridStep) {
            for (double z = -gridRange; z <= gridRange; z += gridStep) {
                Point3D p = {x, y, z};

                // Calculate magnetic field contributions from both coils
                Vector3D B1 = calculateMagneticField(p, coil1, currentCoil1);
                Vector3D B2 = calculateMagneticField(p, coil2, currentCoil2);

                // Total magnetic field is the vector sum of contributions
                Vector3D B = B1 + B2;

                fieldData.push_back({p, B});
            }
        }
    }

    // Save the magnetic field data to a file
    saveFieldToFile("a_coil_pair_magnetic_field_data.txt", fieldData);

    return 0;
}

