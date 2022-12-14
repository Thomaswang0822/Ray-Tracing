#include <glm/glm.hpp>
#include <GLUT/glut.h>
#include <OpenGL/gl3.h>
#include <OpenGL/glext.h>

#include "Ray.h"
#include "Intersection.h"
#include "Camera.h"
#include "Image.h"
#include "RTScene.h"

#ifndef __RAYTRACER_H__
#define __RAYTRACER_H__

namespace RayTracer {
    void Raytrace(Camera* cam, RTScene* scene, Image &image); //page 9
    Ray RayThruPixel(Camera* cam, int i, int j, int width, int height);//page 10,18
    Intersection Intersect(Ray ray, Triangle &triangle); //page 30, 33
    Intersection Intersect(Ray ray, RTScene* scene); //page 11, 28, 31
    glm::vec4 FindColor(Intersection hit, RTScene* scene, int recursion_depth); //page 15

    // global-illumination version
    void RaytraceGlobal(Camera* cam, RTScene* scene, Image &image, int N);
    Ray RandomRayThruPixel(Camera* cam, int i, int j, int width, int height, float p, float q);
    // color by a single random ray: it will be called N times for 1 pixel.
    // Use Russian Roulette
    glm::vec4 FindColorRR(Ray ray, RTScene* scene);
};

#endif