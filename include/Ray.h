#include <glm/glm.hpp>
#include <GLUT/glut.h>
#include <OpenGL/gl3.h>
#include <OpenGL/glext.h>

#ifndef __RAY_H__
#define __RAY_H__
class Ray {
public:
    glm::vec3 p0; // basepoint
    glm::vec3 dir; // direction

    // Member Initialization:
    Ray(glm::vec3 basepoint, glm::vec3 direction)
        : p0(basepoint), dir(direction) {}
    // But can still use the default (trivial) Constructor
    Ray() = default;
};
#endif