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

int main (){
     // Parameters for Maxwell's coil pair
    double majorRadius = 1.5;    // Major radius of the torus
    double minorRadius = 0.1;    // Minor radius (thickness of the torus)
    double separation = 3.0;     // Distance between the two torus coils
    int numMajorPoints = 500;    // Points along the major circle
    int numMinorPoints = 500;     // Points along the minor circle

    // Generate points for the first torus
    std::vector<Point3D> coil1 = generateTorus(majorRadius, minorRadius, -separation / 2.0, 0.0, 0.0, numMajorPoints, numMinorPoints);

    // Generate points for the second torus
    std::vector<Point3D> coil2 = generateTorus(majorRadius, minorRadius, separation / 2.0, 0.0, 0.0, numMajorPoints, numMinorPoints);

    // Generate poins for the third torus (-z-axis)
    std::vector<Point3D> coil3 = generateTorus_zaxis(majorRadius, minorRadius, 0.0, 0.0, -2.0, numMajorPoints, numMinorPoints);

    // Generate points for the fourth torus (z-axis)
    std::vector<Point3D> coil4 = generateTorus_zaxis(majorRadius, minorRadius, 0.0, 0.0, 2.0, numMajorPoints, numMinorPoints);

    // Save coil geometry
    saveCoilGeometryToFile("Maxwell_coil_model.txt", coil1, coil2, coil3, coil4);

    return 0;
}