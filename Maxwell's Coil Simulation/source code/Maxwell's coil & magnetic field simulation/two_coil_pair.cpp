#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <tuple>
#include <omp.h>
#include <iomanip> // For setting output precision

// Structure to hold 3D point data
struct Point3D {
    double x, y, z;
};

// Structure to hold 3D vector data
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

// Function to generate points for a torus-shaped coil
std::vector<Point3D> generateTorus(double majorRadius, double minorRadius, double centerX, double centerY, double centerZ, int numMajorPoints, int numMinorPoints) {
    std::vector<Point3D> coilPoints;
    double majorAngleStep = 2 * M_PI / numMajorPoints;
    double minorAngleStep = 2 * M_PI / numMinorPoints;

    for (int i = 0; i < numMajorPoints; ++i) {
        double majorAngle = i * majorAngleStep;

        // Center of the minor circle in the torus
        double baseX = centerX;
        double baseY = centerY + majorRadius * cos(majorAngle);
        double baseZ = centerZ + majorRadius * sin(majorAngle);

        for (int j = 0; j < numMinorPoints; ++j) {
            double minorAngle = j * minorAngleStep;

            // Points distributed around the minor circle
            double offsetX = minorRadius * cos(minorAngle);
            double offsetY = minorRadius * sin(minorAngle) * cos(majorAngle);
            double offsetZ = minorRadius * sin(minorAngle) * sin(majorAngle);

            coilPoints.push_back({baseX + offsetX, baseY + offsetY, baseZ + offsetZ});
        }
    }

    return coilPoints;
}

std::vector<Point3D> generateTorus_zaxis(double majorRadius, double minorRadius, double centerX, double centerY, double centerZ, int numMajorPoints, int numMinorPoints) {
    std::vector<Point3D> coilPoints;
    double majorAngleStep = 2 * M_PI / numMajorPoints;
    double minorAngleStep = 2 * M_PI / numMinorPoints;

    for (int i = 0; i < numMajorPoints; ++i) {
        double majorAngle = i * majorAngleStep;

        // Center of the minor circle in the torus
        double baseX = centerX + majorRadius * sin(majorAngle);
        double baseY = centerY + majorRadius * cos(majorAngle);
        double baseZ = centerZ;

        for (int j = 0; j < numMinorPoints; ++j) {
            double minorAngle = j * minorAngleStep;

            // Points distributed around the minor circle
            double offsetX = minorRadius * sin(minorAngle) * sin(majorAngle);
            double offsetY = minorRadius * sin(minorAngle) * cos(majorAngle);
            double offsetZ = minorRadius * cos(minorAngle);

            coilPoints.push_back({baseX + offsetX, baseY + offsetY, baseZ + offsetZ});
        }
    }

    return coilPoints;
}

// Function to save coil geometry to a file
void saveCoilGeometryToFile(const std::string& filename, const std::vector<Point3D>& coil1, const std::vector<Point3D>& coil2, 
                            const std::vector<Point3D>& coil3, const std::vector<Point3D>& coil4) {
    std::ofstream outFile(filename);

    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    outFile << std::setprecision(10); // Set precision for output

    outFile << "Coil 1 Points:\n";
    for (const auto& point : coil1) {
        outFile << point.x << ", " << point.y << ", " << point.z << "\n";
    }

    outFile << "\nCoil 2 Points:\n";
    for (const auto& point : coil2) {
        outFile << point.x << ", " << point.y << ", " << point.z << "\n";
    }

    outFile << "\nCoil 3 Points:\n";
    for (const auto& point : coil3) {
        outFile << point.x << ", " << point.y << ", " << point.z << "\n";
    }

    outFile << "\nCoil 4 Points:\n";
    for (const auto& point : coil4) {
        outFile << point.x << ", " << point.y << ", " << point.z << "\n";
    }

    outFile.close();
    std::cout << "Coil geometry saved to " << filename << std::endl;
}

// Function to calculate the magnetic field at a point using Biot-Savart's Law
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

        double r_mag2 = rx * rx + ry * ry + rz * rz;
        double r_mag = std::sqrt(r_mag2);
        if (r_mag < 1e-9) continue; // Avoid division by zero or near-zero values

        double r_mag3 = r_mag2 * r_mag;

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

    outFile << std::setprecision(10); // Set precision for output

    for (const auto& data : fieldData) {
        outFile << data.first.x << ", " << data.first.y << ", " << data.first.z << ", "
                << data.second.x << ", " << data.second.y << ", " << data.second.z << "\n";
    }

    outFile.close();
    std::cout << "Magnetic field data saved to " << filename << std::endl;
}

int main (){
     // Parameters for Maxwell's coil pair
    double majorRadius = 3.0;    // Major radius of the torus
    double minorRadius = 0.1;    // Minor radius (thickness of the torus)
    double separation = 7.0;     // Distance between the two torus coils
    int numMajorPoints = 200;    // Points along the major circle
    int numMinorPoints = 200;     // Points along the minor circle
    double current = 1.0;

    // Generate points for the first torus
    std::vector<Point3D> coil1 = generateTorus(majorRadius, minorRadius, separation / 2.0, 0.0, 0.0, numMajorPoints, numMinorPoints);

    // Generate points for the second torus
    std::vector<Point3D> coil2 = generateTorus(majorRadius, minorRadius, -separation / 2.0, 0.0, 0.0, numMajorPoints, numMinorPoints);

    // Generate points for the third torus (-z-axis)
    std::vector<Point3D> coil3 = generateTorus_zaxis(majorRadius, minorRadius, 0.0, 0.0, separation / 2.0, numMajorPoints, numMinorPoints);

    // Generate points for the fourth torus (z-axis)
    std::vector<Point3D> coil4 = generateTorus_zaxis(majorRadius, minorRadius, 0.0, 0.0, -separation / 2.0, numMajorPoints, numMinorPoints);

    // Combine the points from all coils
    std::vector<Point3D> allCoils;
    allCoils.insert(allCoils.end(), coil1.begin(), coil1.end());
    allCoils.insert(allCoils.end(), coil2.begin(), coil2.end());
    allCoils.insert(allCoils.end(), coil3.begin(), coil3.end());
    allCoils.insert(allCoils.end(), coil4.begin(), coil4.end());

    // Save coil geometry
    saveCoilGeometryToFile("4_maxwell_coil_model.txt", coil1, coil2, coil3, coil4);

    // Define a grid of points to calculate the magnetic field
    std::vector<std::pair<Point3D, Vector3D>> fieldData;
    double gridRange = 7.0; // Range of the grid (from -gridRange to +gridRange)
    double gridStep = 1.0; // Step size for the grid

    // Parallel computation of the magnetic field
    #pragma omp parallel
    {
        std::vector<std::pair<Point3D, Vector3D>> localFieldData;

        #pragma omp for collapse(3) schedule(dynamic)
        for (int ix = 0; ix <= static_cast<int>(2 * gridRange / gridStep); ++ix) {
            for (int iy = 0; iy <= static_cast<int>(2 * gridRange / gridStep); ++iy) {
                for (int iz = 0; iz <= static_cast<int>(2 * gridRange / gridStep); ++iz) {
                    double x = -gridRange + ix * gridStep;
                    double y = -gridRange + iy * gridStep;
                    double z = -gridRange + iz * gridStep;

                    Point3D p = {x, y, z};
                    Vector3D B = calculateMagneticField(p, allCoils, current);
                    localFieldData.push_back({p, B});
                }
            }
        }

        #pragma omp critical
        fieldData.insert(fieldData.end(), localFieldData.begin(), localFieldData.end());
    }

    // Save the magnetic field data to a file
    saveFieldToFile("4_coil_magnetic_field_data.txt", fieldData);

    return 0;
}
