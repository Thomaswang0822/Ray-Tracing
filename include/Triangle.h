#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <vector>
#include "Material.h"
#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__
struct Triangle {
    std::vector<glm::vec3> P; // 3 positions
    std::vector<glm::vec3> N; // 3 normals
    Material* material = NULL;

    /* This function is called in Scene.cpp DFS loop;
     * a copy of Triangle from that in Geometry::elements will call this;
     * so we need to UPDATE values in P and N 
     */
    void transPN(glm::mat4 tr) {  // given transformation matrix
        for (auto &pos : P) {
            pos = glm::vec3( tr * glm::vec4(pos, 1.0f) );
        }
        glm::mat4 invtr = glm::inverse(tr); // used for normal vector
        for (auto &normal : N) {
            normal = glm::mat3(invtr) * normal;
            normal = glm::normalize(normal);    // unnecessary in theory; just in case
        }
    }
};
#endif