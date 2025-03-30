#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    float x, y, z;
} Vertex;

// Struktur Face diperluas dengan atribut 'color': 
// 0 = pola chess board (misal: darkgreen), 1 = pola chess board (misal: lightgreen),
// 2 = ribbon (misalnya, akan di-plot dengan warna khusus, misal orange)
typedef struct {
    int v1, v2, v3;
    int color;
} Face;

// Fungsi transformasi: skala pada sumbu Y dan rotasi 45° terhadap sumbu Y
Vertex transformVertex(Vertex v, float scaleY, float cosTheta, float sinTheta) {
    Vertex vt;
    vt.y = v.y * scaleY;
    vt.x = v.x * cosTheta + v.z * sinTheta;
    vt.z = -v.x * sinTheta + v.z * cosTheta;
    return vt;
}

int main() {
    // --- GEOMETRI KUBUS DENGAN CHESS BOARD ---
    int subdiv = 4; // subdivisi per sisi (grid 4x4) untuk kubus
    int verticesPerFace = (subdiv + 1) * (subdiv + 1);
    int totalCubeVertices = 6 * verticesPerFace;
    int facesPerFace = subdiv * subdiv * 2; // 2 segitiga per sel
    int totalCubeFaces = 6 * facesPerFace;
    
    // --- GEOMETRI EXTRA: 3 pita (ribbon)
    // Satu pita tetap pada front face (ribbon 0)
    // Dua pita pada face kiri (sumbu x negatif): ribbon 1 (atas) dan ribbon 2 (bawah)
    int ribbonCount = 4;
    int verticesPerRibbon = 4; // tiap ribbon berupa rectangle dengan 4 vertex
    int facesPerRibbon = 2;    // tiap ribbon dibagi menjadi 2 segitiga
    int extraVerticesCount = ribbonCount * verticesPerRibbon; // 3 * 4 = 12
    int extraFacesCount = ribbonCount * facesPerRibbon;           // 3 * 2 = 6

    // Total geometry baru:
    int totalVertices = totalCubeVertices + extraVerticesCount;
    int totalFaces = totalCubeFaces + extraFacesCount;
    
    // Alokasikan array untuk menyimpan vertex dan face
    Vertex *vertices = (Vertex *)malloc(totalVertices * sizeof(Vertex));
    Face *faces = (Face *)malloc(totalFaces * sizeof(Face));
    if (!vertices || !faces) {
        printf("Gagal mengalokasikan memori!\n");
        return 1;
    }
    
    float step = 2.0f / subdiv; // tiap langkah dari -1 ke 1
    int vertexIndex = 0;
    
    // Generate vertex untuk masing-masing sisi kubus
    // Urutan sisi: 0-front, 1-back, 2-left, 3-right, 4-top, 5-bottom.
    for (int f = 0; f < 6; f++) {
        for (int i = 0; i <= subdiv; i++) {
            float v = -1.0f + i * step;
            for (int j = 0; j <= subdiv; j++) {
                float u = -1.0f + j * step;
                Vertex vert;
                switch (f) {
                    case 0: // Front: z = 1
                        vert.x = u;
                        vert.y = v;
                        vert.z = 1.0f;
                        break;
                    case 1: // Back: z = -1, balik nilai u
                        vert.x = -u;
                        vert.y = v;
                        vert.z = -1.0f;
                        break;
                    case 2: // Left: x = -1
                        vert.x = -1.0f;
                        vert.y = v;
                        vert.z = u;
                        break;
                    case 3: // Right: x = 1, balik u
                        vert.x = 1.0f;
                        vert.y = v;
                        vert.z = -u;
                        break;
                    case 4: // Top: y = 1, gunakan u untuk x dan v untuk -z
                        vert.x = u;
                        vert.y = 1.0f;
                        vert.z = -v;
                        break;
                    case 5: // Bottom: y = -1
                        vert.x = u;
                        vert.y = -1.0f;
                        vert.z = v;
                        break;
                    default:
                        vert.x = vert.y = vert.z = 0.0f;
                        break;
                }
                vertices[vertexIndex++] = vert;
            }
        }
    }
    
    // Parameter transformasi
    float scaleY = 0.5f;  // mempersempit sumbu Y
    float angle = 45.0f * (3.14159265f / 180.0f);
    float cosTheta = cosf(angle);
    float sinTheta = sinf(angle);
    
    // Terapkan transformasi ke semua vertex kubus
    for (int i = 0; i < totalCubeVertices; i++) {
        vertices[i] = transformVertex(vertices[i], scaleY, cosTheta, sinTheta);
    }
    
    // Generate face untuk kubus (polanya chess board)
    int faceIndex = 0;
    int faceVertexOffset = 0;
    for (int f = 0; f < 6; f++) {
        for (int i = 0; i < subdiv; i++) {
            for (int j = 0; j < subdiv; j++) {
                int idx0 = faceVertexOffset + i * (subdiv + 1) + j;
                int idx1 = faceVertexOffset + i * (subdiv + 1) + (j + 1);
                int idx2 = faceVertexOffset + (i + 1) * (subdiv + 1) + (j + 1);
                int idx3 = faceVertexOffset + (i + 1) * (subdiv + 1) + j;
                int color = ((i + j) % 2 == 0) ? 0 : 1;
                faces[faceIndex].v1 = idx0;
                faces[faceIndex].v2 = idx1;
                faces[faceIndex].v3 = idx2;
                faces[faceIndex].color = color;
                faceIndex++;
                faces[faceIndex].v1 = idx0;
                faces[faceIndex].v2 = idx2;
                faces[faceIndex].v3 = idx3;
                faces[faceIndex].color = color;
                faceIndex++;
            }
        }
        faceVertexOffset += verticesPerFace;
    }
    
    // --- Generate Extra Geometry: 3 Ribbons ---
    // Kita buat:
    // Ribbon 0 (Front Ribbon): terletak pada front face, dengan tepi kiri pada x = -1
    //   - Base rectangle (sistem asli): 
    //       v0: (-1, -0.4, 1), v1: (-0.5, -0.4, 1), v2: (-0.5, 0.4, 1), v3: (-1, 0.4, 1)
    //   - Ekstrusi: tambahkan thickness pada sumbu z (ke depan)
    // Ribbon 1 (Top Left Ribbon): terletak pada left face (x = -1), di bagian atas
    //   - Base rectangle (sistem asli):
    //       v0: (-1, 0.0, 1), v1: (-1, 0.0, 1)  (sudut yang sama untuk sisi yang bersinggungan dengan front ribbon)
    //       v2: (-1, 0.0, 1) tidak, kita atur agar berbentuk rectangle dengan dua kolom
    //   - Kita definisikan sebagai:
    //       Kolom 0: berada di x = -1, dengan y dari 0.0 ke 0.4, z = 1
    //       Kolom 1: berada di x = -1 - thickness, dengan y sama, z = 1
    // Ribbon 2 (Bottom Left Ribbon): terletak pada left face (x = -1), di bagian bawah
    //       Kolom 0: x = -1, dengan y dari -0.4 ke 0.0, z = 1
    //       Kolom 1: x = -1 - thickness, dengan y sama, z = 1
    float thickness = 0.2f; // ekstrusi (untuk front ribbon: ke sumbu z; untuk left ribbons: ke negatif x)
    
    int extraVertexStart = totalCubeVertices; // indeks awal extra vertex
    // Ribbon 0: Front Ribbon
    {
        int baseIdx = extraVertexStart + 0; // ribbon 2 dimulai dari sini
        // Gunakan y dari -0.4 ke 0.0, x sama seperti ribbon 1, z = 1
        Vertex col0[2] = {
            { -1.0f, 0.0f, 1.0f },
            { -1.0f,  0.4f, 1.0f }
        };
        Vertex col1[2] = {
            { -1.0f, 0.0f, 2.0f + thickness },
            { -1.0f,  0.4f, 2.0f + thickness }
        };
        vertices[baseIdx + 0] = transformVertex(col0[0], scaleY, cosTheta, sinTheta);
        vertices[baseIdx + 1] = transformVertex(col1[0], scaleY, cosTheta, sinTheta);
        vertices[baseIdx + 2] = transformVertex(col1[1], scaleY, cosTheta, sinTheta);
        vertices[baseIdx + 3] = transformVertex(col0[1], scaleY, cosTheta, sinTheta);
        
        // Buat 2 face untuk ribbon 2
        faces[faceIndex].v1 = baseIdx + 0;
        faces[faceIndex].v2 = baseIdx + 1;
        faces[faceIndex].v3 = baseIdx + 2;
        faces[faceIndex].color = 2;
        faceIndex++;
        faces[faceIndex].v1 = baseIdx + 0;
        faces[faceIndex].v2 = baseIdx + 2;
        faces[faceIndex].v3 = baseIdx + 3;
        faces[faceIndex].color = 2;
        faceIndex++;
    }
    
    // Ribbon 1: Top Left Ribbon (di left face, bagian atas)
    {
        int baseIdx = extraVertexStart + 4; // ribbon 1 dimulai dari indeks ini
        // Pada left face, x = -1 untuk kolom 0, x = -1 - thickness untuk kolom 1
        // Gunakan y dari 0.0 ke 0.4, dan z = 1 (sama dengan front)
        Vertex col0[2] = {
            { -1.0f,  -0.4f, 1.0f },
            { -1.0f,  0.0f, 1.0f }
        };
        Vertex col1[2] = {
            { -1.0f, -0.4f, 2.0f + thickness},
            { -1.0f, 0.0f, 2.0f + thickness }
        };
        // Simpan ke array extra vertex untuk ribbon 1 dengan urutan: 
        // baris bawah: col0[0], col1[0]; baris atas: col0[1], col1[1]
        vertices[baseIdx + 0] = transformVertex(col0[0], scaleY, cosTheta, sinTheta);
        vertices[baseIdx + 1] = transformVertex(col1[0], scaleY, cosTheta, sinTheta);
        vertices[baseIdx + 2] = transformVertex(col1[1], scaleY, cosTheta, sinTheta);
        vertices[baseIdx + 3] = transformVertex(col0[1], scaleY, cosTheta, sinTheta);
        
        // Buat 2 face untuk ribbon 1
        faces[faceIndex].v1 = baseIdx + 0;
        faces[faceIndex].v2 = baseIdx + 1;
        faces[faceIndex].v3 = baseIdx + 2;
        faces[faceIndex].color = 2;
        faceIndex++;
        faces[faceIndex].v1 = baseIdx + 0;
        faces[faceIndex].v2 = baseIdx + 2;
        faces[faceIndex].v3 = baseIdx + 3;
        faces[faceIndex].color = 2;
        faceIndex++;
    }
    
    // Ribbon 2: Bottom Left Ribbon (di left face, bagian bawah)
    {
        int baseIdx = extraVertexStart + 8; // ribbon 2 dimulai dari sini
        // Gunakan y dari -0.4 ke 0.0, x sama seperti ribbon 1, z = 1
        Vertex col0[2] = {
            { -1.0f, 0.0f, 1.0f },
            { -1.0f,  0.4f, 1.0f }
        };
        Vertex col1[2] = {
            { -2.0f - thickness, 0.0f, 1.0f },
            { -2.0f - thickness,  0.4f, 1.0f }
        };
        vertices[baseIdx + 0] = transformVertex(col0[0], scaleY, cosTheta, sinTheta);
        vertices[baseIdx + 1] = transformVertex(col1[0], scaleY, cosTheta, sinTheta);
        vertices[baseIdx + 2] = transformVertex(col1[1], scaleY, cosTheta, sinTheta);
        vertices[baseIdx + 3] = transformVertex(col0[1], scaleY, cosTheta, sinTheta);
        
        // Buat 2 face untuk ribbon 2
        faces[faceIndex].v1 = baseIdx + 0;
        faces[faceIndex].v2 = baseIdx + 1;
        faces[faceIndex].v3 = baseIdx + 2;
        faces[faceIndex].color = 2;
        faceIndex++;
        faces[faceIndex].v1 = baseIdx + 0;
        faces[faceIndex].v2 = baseIdx + 2;
        faces[faceIndex].v3 = baseIdx + 3;
        faces[faceIndex].color = 2;
        faceIndex++;
    }
    {
        int baseIdx = extraVertexStart + 12; // ribbon 2 dimulai dari sini
        // Gunakan y dari -0.4 ke 0.0, x sama seperti ribbon 1, z = 1
        Vertex col0[2] = {
            { -1.0f, -0.4f, 1.0f },
            { -1.0f, 0.0f, 1.0f }
        };
        Vertex col1[2] = {
            { -2.0f - thickness, -0.4f, 1.0f },
            { -2.0f - thickness, 0.0f, 1.0f }
        };
        vertices[baseIdx + 0] = transformVertex(col0[0], scaleY, cosTheta, sinTheta);
        vertices[baseIdx + 1] = transformVertex(col1[0], scaleY, cosTheta, sinTheta);
        vertices[baseIdx + 2] = transformVertex(col1[1], scaleY, cosTheta, sinTheta);
        vertices[baseIdx + 3] = transformVertex(col0[1], scaleY, cosTheta, sinTheta);
        
        // Buat 2 face untuk ribbon 2
        faces[faceIndex].v1 = baseIdx + 0;
        faces[faceIndex].v2 = baseIdx + 1;
        faces[faceIndex].v3 = baseIdx + 2;
        faces[faceIndex].color = 2;
        faceIndex++;
        faces[faceIndex].v1 = baseIdx + 0;
        faces[faceIndex].v2 = baseIdx + 2;
        faces[faceIndex].v3 = baseIdx + 3;
        faces[faceIndex].color = 2;
        faceIndex++;
    }
    // Simpan data model ke file biner
    FILE *fp = fopen("cube_chess_model.bin", "wb");
    if (fp == NULL) {
        printf("Gagal membuka file untuk menulis!\n");
        free(vertices);
        free(faces);
        return 1;
    }
    
    fwrite(&totalVertices, sizeof(int), 1, fp);
    fwrite(vertices, sizeof(Vertex), totalVertices, fp);
    fwrite(&totalFaces, sizeof(int), 1, fp);
    fwrite(faces, sizeof(Face), totalFaces, fp);
    
    fclose(fp);
    printf("Model kubus dengan chess board, rotasi 45°, lebar Y dipersempit, dan 3 ribbon (1 di sumbu z, 2 di sumbu x negatif yang bersinggungan) telah disimpan dalam file cube_chess_model.bin\n");
    
    free(vertices);
    free(faces);
    return 0;
}
