#include <glm/glm.hpp>
#include <GLUT/glut.h>
#include <OpenGL/gl3.h>
#include <OpenGL/glext.h>
#include "Triangle.h"

#ifndef __INTERSECTION_H__
#define __INTERSECTION_H__
class Intersection {
public:
    glm::vec3 P; // position of the intersection
    glm::vec3 N; // surface normal
    glm::vec3 V; // direction to incoming ray
    Triangle* triangle;
    float dist = INFINITY; // distance to the source of ray

    // Member Initialization:
    Intersection(glm::vec3 p, glm::vec3 n, glm::vec3 v, Triangle* triPt, float d):
        P(p), N(n), V(v), triangle(triPt), dist(d) {};
    // But can still use the default (trivial) Constructor
    // default Intersection is "no-intersection", by dist==INFINITY
    Intersection() = default;
};
#endif