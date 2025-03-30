#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

// Updated Vertex structure: storing Cartesian coordinates and spherical angles.
struct Vertex {
    float x, y, z;
    float theta, phi; // Spherical coordinates (polar and azimuthal angles)
};

// Structure to hold the mesh data: vertices and triangle indices.
struct Mesh {
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
};

// Create a sphere mesh with the given radius, stacks (latitudinal subdivisions),
// and slices (longitudinal subdivisions). Also store the spherical angles.
Mesh createSphereMesh(float radius, unsigned int stacks, unsigned int slices) {
    Mesh mesh;
    
    // Generate vertices using spherical coordinates.
    for (unsigned int i = 0; i <= stacks; i++) {
        float theta = i * M_PI / stacks;
        float sinTheta = std::sin(theta);
        float cosTheta = std::cos(theta);
        
        for (unsigned int j = 0; j <= slices; j++) {
            float phi = j * 2 * M_PI / slices;
            Vertex vertex;
            vertex.theta = theta;
            vertex.phi = phi;
            // Compute Cartesian coordinates for the base sphere.
            vertex.x = radius * sinTheta * std::cos(phi);
            vertex.y = radius * sinTheta * std::sin(phi);
            vertex.z = radius * cosTheta;
            mesh.vertices.push_back(vertex);
        }
    }
    
    // Generate triangle indices (each rectangular patch is divided into two triangles).
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

// Propagating pulse with reflection and damping: computes a traveling Gaussian pulse
// that reflects at the antipodal point and whose amplitude decays over time.
// gamma is the angular distance between the vertex and the source.
// The pulse displacement is computed as:
//   Δr = A * exp( - d^2/(2σ^2) ) * cos(2π f d) * exp(-decayRate * t)
// where d = γ - d_eff, and d_eff is the effective travel distance (triangle-wave behavior).
float propagatingPulseReflecting(float theta, float phi, float time, float amplitude,
                                 float thetaSource, float phiSource, float sigma,
                                 float waveSpeed, float frequency, float decayRate) {
    // Compute angular distance gamma using the spherical law of cosines.
    float cosGamma = std::sin(theta) * std::sin(thetaSource) * std::cos(phi - phiSource)
                   + std::cos(theta) * std::cos(thetaSource);
    if(cosGamma > 1.0f) cosGamma = 1.0f;
    if(cosGamma < -1.0f) cosGamma = -1.0f;
    float gamma = std::acos(cosGamma);
    
    // Compute the travel distance with reflection.
    float totalTravel = waveSpeed * time;
    float modTravel = std::fmod(totalTravel, 2 * M_PI);
    float effectiveTravel = (modTravel > M_PI) ? (2 * M_PI - modTravel) : modTravel;
    
    // Compute difference between vertex angular distance and effective pulse front.
    float d = gamma - effectiveTravel;
    
    // Gaussian envelope and oscillatory factor.
    float envelope = std::exp( - (d * d) / (2 * sigma * sigma) );
    float oscillation = std::cos(2 * M_PI * frequency * d);
    
    // Apply damping so that when the wave meets itself the amplitude decays.
    float dampingFactor = std::exp(-decayRate * time);
    
    return amplitude * dampingFactor * envelope * oscillation;
}

// Update vertex positions by summing the displacement from the first pulse
// and, when time >= halfTime, the displacement from a second pulse.
// The second pulse is generated at the same source location as the first pulse,
// but its time is offset by halfTime.
void updateVertexPositions(Mesh &mesh, float baseRadius,
                           // First pulse parameters:
                           float amplitude1, float time,
                           float thetaSource1, float phiSource1, float sigma1,
                           float waveSpeed1, float frequency1, float decayRate1,
                           // Second pulse parameters (same source location):
                           float amplitude2, float sigma2,
                           float waveSpeed2, float frequency2, float decayRate2,
                           float halfTime) {
    for (auto &v : mesh.vertices) {
        float disp1 = propagatingPulseReflecting(v.theta, v.phi, time, amplitude1,
                                                 thetaSource1, phiSource1, sigma1,
                                                 waveSpeed1, frequency1, decayRate1);
        float disp2 = 0.0f;
        if (time >= halfTime) {
            // For the second pulse, subtract halfTime from current time so that it starts at t = 0.
            float t2 = time - halfTime;
            disp2 = propagatingPulseReflecting(v.theta, v.phi, t2, amplitude2,
                                               thetaSource1, phiSource1, sigma2,
                                               waveSpeed2, frequency2, decayRate2);
        }
        float newRadius = baseRadius + disp1 + disp2;
        // Update Cartesian coordinates.
        v.x = newRadius * std::sin(v.theta) * std::cos(v.phi);
        v.y = newRadius * std::sin(v.theta) * std::sin(v.phi);
        v.z = newRadius * std::cos(v.theta);
    }
}

int main() {
    // Parameters for the base sphere.
    float baseRadius = 5.0f;
    unsigned int stacks = 50;       // Latitudinal subdivisions.
    unsigned int slices = 50;       // Longitudinal subdivisions.
    
    // Time parameters.
    float timeStart = 0.0f;
    float timeEnd = 6 * M_PI;       // Total simulation time.
    unsigned int numFrames = 1500;   // Total number of frames.
    float dt = (timeEnd - timeStart) / numFrames;
    float halfTime = timeEnd / 4.0f;
    
    // First pulse parameters.
    float pulseAmplitude1 = 2.0f;       // Amplitude of the first pulse.
    float thetaSource1 = M_PI / 2.0f;     // Source located at the equator.
    float phiSource1 = 0.0f;              // Azimuthal position of the first source.
    float sigma1 = 0.2f;                // Spatial width of the first pulse.
    float waveSpeed1 = 0.6f;            // Propagation speed of the first pulse.
    float frequency1 = 2.0f;            // Frequency of oscillation within the first pulse.
    float decayRate1 = 0.15f;            // Decay rate for the first pulse.
    
    // Second pulse parameters (same source as the first pulse).
    float pulseAmplitude2 = 1.5f;       // Amplitude of the second pulse.
    float sigma2 = 0.2f;                // Spatial width of the second pulse.
    float waveSpeed2 = 0.6f;            // Propagation speed of the second pulse.
    float frequency2 = 2.0f;            // Frequency of oscillation within the second pulse.
    float decayRate2 = 0.3f;            // Decay rate for the second pulse.
    
    // Create the base sphere mesh.
    Mesh sphereMesh = createSphereMesh(baseRadius, stacks, slices);
    
    // Open a file to write JSON animation data.
    std::ofstream outFile("animation_data.json");
    if (!outFile) {
        std::cerr << "Error opening file for writing." << std::endl;
        return 1;
    }
    
    // Write JSON header with fixed floating-point format.
    outFile << std::fixed << std::setprecision(6);
    outFile << "{\n";
    
    // Write constant indices (mesh connectivity).
    outFile << "  \"indices\": [\n    ";
    for (size_t i = 0; i < sphereMesh.indices.size(); i++) {
        outFile << sphereMesh.indices[i];
        if (i < sphereMesh.indices.size() - 1)
            outFile << ", ";
        if ((i + 1) % 10 == 0)
            outFile << "\n    ";
    }
    outFile << "\n  ],\n";
    
    // Write frame data.
    outFile << "  \"frames\": [\n";
    for (unsigned int frame = 0; frame <= numFrames; frame++) {
        float currentTime = timeStart + frame * dt;
        // Update vertex positions with the combined pulses.
        updateVertexPositions(sphereMesh, baseRadius,
                              pulseAmplitude1, currentTime,
                              thetaSource1, phiSource1, sigma1,
                              waveSpeed1, frequency1, decayRate1,
                              pulseAmplitude2, sigma2,
                              waveSpeed2, frequency2, decayRate2,
                              halfTime);
        
        outFile << "    {\n";
        outFile << "      \"time\": " << currentTime << ",\n";
        outFile << "      \"vertices\": [\n        ";
        for (size_t i = 0; i < sphereMesh.vertices.size(); i++) {
            const auto &v = sphereMesh.vertices[i];
            outFile << "[" << v.x << ", " << v.y << ", " << v.z << "]";
            if (i < sphereMesh.vertices.size() - 1)
                outFile << ", ";
            if ((i + 1) % 5 == 0)
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
