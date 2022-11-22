/*
Geometry used for Ray Tracing.

Replace traditional rasterization parameter (vertex buffer, index buffer, etc)
with a list of triangles.
*/
#include <vector>
#include "Triangle.h"
#ifndef __RTGEOMETRY_H__
#define __RTGEOMETRY_H__
class RTGeometry {
public:
    int count; // number of elements to draw
    std::vector<Triangle> elements; // list of triangles
    virtual void init(){};
    virtual void init(const char* s){};
};
#endif