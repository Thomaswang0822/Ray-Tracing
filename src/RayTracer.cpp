#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#ifdef __APPLE__
#include <OpenGL/gl3.h>
#include <OpenGL/glext.h>
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

#include "RayTracer.h"
static int max_depth = 5;

void RayTracer::Raytrace(Camera cam, RTScene* scene, Image &image) {
    std::cout << "Raytrace() called" << std::endl;
    int w = image.width; int h = image.height;
    Ray ray;
    Intersection hit; 
    for (int j=0; j<h; j++){
        for (int i=0; i<w; i++){
            ray = RayThruPixel( cam, i, j, w, h );
            // std::cout << "ray found" << std::endl;
            hit = Intersect( ray, scene );
            // std::cout << "hit found" << std::endl;
            // pixels is 1d vector
            image.pixels[j*w + i] = FindColor( hit, max_depth );
            // std::cout << "pixel has color" << std::endl;
            //break; // debug use: compute a single pixel only
        }
    }
}; // page 9

Ray RayTracer::RayThruPixel(Camera cam, int i, int j, int width, int height){
    float alpha = (float(i) + 0.5f) * 2.f/width - 1.f;
    float beta = 1.f - (float(j) + 0.5f) * 2.f/height;
    float fovy_h = 0.5f * cam.fovy / 180.f * M_PI;     // half fovy, turn degree into rad

    glm::vec3 w = glm::normalize(cam.eye - cam.target);
    glm::vec3 u = glm::normalize(glm::cross(cam.up, w));
    glm::vec3 v = glm::cross(w, u);
    glm::vec3 dir = glm::normalize(
        alpha * cam.aspect * tanf(fovy_h) * u 
        + beta * tanf(fovy_h)
        - w
    ); 
    return Ray(cam.eye, dir);   // return Ray in World coord

};//page 10,18

Intersection RayTracer::Intersect(Ray ray, Triangle triangle){
    // extract variable for readbility
    glm::vec3 p1, p2, p3, n1, n2, n3;
    p1 = triangle.P[0]; n1 = triangle.N[0];
    p2 = triangle.P[1]; n2 = triangle.N[1];
    p3 = triangle.P[2]; n3 = triangle.N[2];

    // construct Ax=b to solve for x = [ lambda[1:3], t]
    glm::mat4 A;    // glm mat is column-major
    A[0] = glm::vec4(p1, 1.f);
    A[1] = glm::vec4(p2, 1.f);
    A[2] = glm::vec4(p3, 1.f);
    A[3] = glm::vec4(-ray.dir, 0.f);
    glm::vec4 b = glm::vec4(ray.p0, 1.f);
    glm::vec4 x = glm::inverse(A) * b;
    // if any of lambda or d negative: no-hit
    for (int i=0; i<4; i++) {
        if (x[i] < 0) {
            return Intersection();
        }
    }
    float lam1, lam2, lam3, d;  // unpack x
    lam1 = x[0]; lam2 = x[1]; lam3 = x[2]; d = x[3];

    glm::vec3 q = lam1*p1 + lam2*p2 +lam3*p3;
    glm::vec3 n = glm::normalize(lam1*n1 + lam2*n2 +lam3*n3);

    return Intersection(q, n, -ray.dir, &triangle, d);
}; //page 30, 33

Intersection RayTracer::Intersect(Ray ray, RTScene* scene){
    float mindist = INFINITY;
    Intersection hit;
    Intersection hit_temp;
    // Find closest intersection; test all triangles
    for (Triangle tri : scene->triangle_soup) {
        hit_temp = Intersect(ray, tri);
        if (hit_temp.dist< mindist){ // closer than previous hit
            mindist = hit_temp.dist;
            hit = hit_temp;
        }
    }
    return hit;
}; //page 11, 28, 31

glm::vec3 RayTracer::FindColor(Intersection hit, int recursion_depth){
    glm::vec3 color = glm::vec3(0.f);
    if (hit.dist < INFINITY) {
        // if hit, white pixel
        color = glm::vec3(1.f);
    }
    // else black
    return color;
}; //page 15
