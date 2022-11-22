/**************************************************
Cube is subclass class of Geometry
that represents a 3D cube.
*****************************************************/
#include "RTGeometry.h"
#include <GLUT/glut.h>
#include <OpenGL/gl3.h>
#include <OpenGL/glext.h>
#ifndef __CUBE_H__
#define __CUBE_H__

class RTCube : public RTGeometry {
public:

    void init(void){
        // positive xyz: right, up, back(further away)
        const GLfloat positions[24][3] ={
            // Front face
            { -0.5f, -0.5f, 0.5f },{ -0.5f, 0.5f, 0.5f },{ 0.5f, 0.5f, 0.5f },{ 0.5f, -0.5f, 0.5f },
            // Back face
            { -0.5f, -0.5f, -0.5f },{ -0.5f, 0.5f, -0.5f },{ 0.5f, 0.5f, -0.5f },{ 0.5f, -0.5f, -0.5f },
            // Left face
            { -0.5f, -0.5f, 0.5f },{ -0.5f, 0.5f, 0.5f },{ -0.5f, 0.5f, -0.5f },{ -0.5f, -0.5f, -0.5f },
            // Right face
            { 0.5f, -0.5f, 0.5f },{ 0.5f, 0.5f, 0.5f },{ 0.5f, 0.5f, -0.5f },{ 0.5f, -0.5f, -0.5f },
            // Top face
            { 0.5f, 0.5f, 0.5f },{ -0.5f, 0.5f, 0.5f },{ -0.5f, 0.5f, -0.5f },{ 0.5f, 0.5f, -0.5f },
            // Bottom face
            { 0.5f, -0.5f, 0.5f },{ -0.5f, -0.5f, 0.5f },{ -0.5f, -0.5f, -0.5f },{ 0.5f, -0.5f, -0.5f }
        };
        // vertex normals
        const GLfloat normals[24][3] = {
            // Front face
            { 0.0f, 0.0f, 1.0f },{ 0.0f, 0.0f, 1.0f },{ 0.0f, 0.0f, 1.0f },{ 0.0f, 0.0f, 1.0f },
            // Back face
            { 0.0f, 0.0f, -1.0f },{ 0.0f, 0.0f, -1.0f },{ 0.0f, 0.0f, -1.0f },{ 0.0f, 0.0f, -1.0f },
            // Left face
            { -1.0f, 0.0f, 0.0f },{ -1.0f, 0.0f, 0.0f },{ -1.0f, 0.0f, 0.0f },{ -1.0f, 0.0f, 0.0f },
            // Right face
            { 1.0f, 0.0f, 0.0f },{ 1.0f, 0.0f, 0.0f },{ 1.0f, 0.0f, 0.0f },{ 1.0f, 0.0f, 0.0f },
            // Top face
            { 0.0f, 1.0f, 0.0f },{ 0.0f, 1.0f, 0.0f },{ 0.0f, 1.0f, 0.0f },{ 0.0f, 1.0f, 0.0f },
            // Bottom face
            { 0.0f, -1.0f, 0.0f },{ 0.0f, -1.0f, 0.0f },{ 0.0f, -1.0f, 0.0f },{ 0.0f, -1.0f, 0.0f }
        };
        // Cube indices
        const GLuint indices[36] = {
            0, 1, 2, 0, 2, 3, // Front face
            4, 5, 6, 4, 6, 7, // Back face
            8, 9, 10, 8, 10, 11, // Left face
            12, 13, 14, 12, 14, 15, // Right face
            16, 17, 18, 16, 18, 19, // Top face
            20, 21, 22, 20, 22, 23 // Bottom face
        };

        // Populate elements: a list of triangle
        // Make use of above positions & normal & indices arrays
        int idx, idx1, idx2;
        for (int i=0; i<36; i+=3) {
            Triangle tri = Triangle();
            // map i to 3 indices
            idx = indices[i]; idx1 = indices[i+1]; idx2 = indices[i+2];
            // each 3 indices is a triangle
            tri.P.push_back( glm::vec3(positions[idx][0], positions[idx][1], positions[idx][2]) );
            tri.P.push_back( glm::vec3(positions[idx1][0], positions[idx1][1], positions[idx1][2]) );
            tri.P.push_back( glm::vec3(positions[idx2][0], positions[idx2][1], positions[idx2][2]) );
            
            tri.N.push_back( glm::vec3(normals[idx][0], normals[idx][1], normals[idx][2]) );
            tri.N.push_back( glm::vec3(normals[idx1][0], normals[idx1][1], normals[idx1][2]) );
            tri.N.push_back( glm::vec3(normals[idx2][0], normals[idx2][1], normals[idx2][2]) );

            elements.push_back(tri);
        }
        
        // count = sizeof(indices)/sizeof(indices[0]);
        count = 6*2;    // 6 faces, 2 triangles each
    }
    
    
};

#endif 
