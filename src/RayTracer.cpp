#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <random>
#include "glm/gtx/string_cast.hpp"
#ifdef __APPLE__
#include <OpenGL/gl3.h>
#include <OpenGL/glext.h>
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

#include "RayTracer.h"
static const int max_depth = 3;

void RayTracer::Raytrace(Camera* cam, RTScene* scene, Image &image) {
    std::cout << "Raytrace() called" << std::endl;
    int w = image.width; int h = image.height;
    Ray ray;
    Intersection hit; 
    // image.pixels.clear();
    #pragma omp parallel for
    for (int j=0; j<h; j++){
        for (int i=0; i<w; i++){
            ray = RayThruPixel( cam, i, j, w, h );
            hit = Intersect( ray, scene );
            // pixels is 1d vector
            image.pixels[(h-1-j)*w + i] = glm::vec3(FindColor(hit, scene, 1));
            // image.pixels.push_back(glm::vec3(FindColor(hit, scene, 1)));
        }
    }
}; // page 9

Ray RayTracer::RayThruPixel(Camera* cam, int i, int j, int width, int height){
    float alpha = (float(i) + 0.5f) * 2.f/width - 1.f;
    float beta = 1.f - (float(j) + 0.5f) * 2.f/height;
    float fovy_h = 0.5f * cam->fovy / 180.f * M_PI;     // half fovy, turn degree into rad

    glm::vec3 w = glm::normalize(cam->eye - cam->target);
    glm::vec3 u = glm::normalize(glm::cross(cam->up, w));
    glm::vec3 v = glm::cross(w, u);
    glm::vec3 dir = glm::normalize(
        alpha * cam->aspect * tanf(fovy_h) * u 
        + beta * tanf(fovy_h) * v
        - w
    ); 
    return Ray(cam->eye, dir);   // return Ray in World coord

};//page 10,18

Intersection RayTracer::Intersect(Ray ray, Triangle& triangle){
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
    for (Triangle &tri : scene->triangle_soup) {
        hit_temp = Intersect(ray, tri);
        if (hit_temp.dist< mindist){ // closer than previous hit
            mindist = hit_temp.dist;
            hit = hit_temp;
        }
    }
    return hit;
}; //page 11, 28, 31

glm::vec4 RayTracer::FindColor(Intersection hit, RTScene* scene, int recursion_depth){
    if (hit.dist == INFINITY) {
        return glm::vec4(0.0f); // if not hit, black background / no light color
    }
    // everything in world coord
    Material* ma = hit.triangle->material;
    glm::vec4 color = ma->emision;
    // glm::vec4 color = glm::vec4(0.0f);
    // For every light
    for (int j=0; j<scene->shader->nlights; j++) {
        // generate 2nd ray
        glm::vec4 lightPos = scene->shader->lightpositions[j];
        // use work-around from HW3, to deal with vec4 light with w = 0 (infinite dist, sunlight)
        glm::vec3 l = glm::normalize(glm::vec3(lightPos) * 1.0f - hit.P * lightPos.w);
        //To avoid self-shadowing, the secondary ray is shot off slightly above the hitting point.
        Ray rayLight = Ray(hit.P+1e-3f*hit.N, l);   // basepoint, dir to light (light - pos)
        // Determine visibility: dist = inf means not visible
        Intersection hitL = Intersect(rayLight, scene);
        int visible = hitL.dist==INFINITY? 0:1;

        // Shading Model; adepted from lighting.frag
        glm::vec4 inside = ma->ambient;
        // Add diffuse part
        inside += ma->diffuse * std::max(glm::dot(hit.N, l), 0.f) * float(visible);

        // Specular part no recursion: 
        // old Blinn–Phong specular reflection if depth reaches limit
        if (recursion_depth == max_depth) {
        // if (recursion_depth<max_depth) {    // DEBUG
            // find h
            glm::vec3 h = glm::normalize(hit.V + l);
            inside += ma->specular * pow( 
                glm::max(glm::dot(hit.N, h), 0.f), 
                ma->shininess
            );
        }

        // Add contribution of this light to color
        color += inside * scene->shader->lightcolors[j];
    }
    // recurvise specular 
    if((recursion_depth < max_depth)) {
        // mirror reflection ray
        Ray ray2 = Ray(hit.P+1e-3f*hit.N, 2.f * glm::dot(hit.N, hit.V)*hit.N - hit.V);    
        Intersection hit2 = Intersect(ray2, scene);
        color += ma->specular * FindColor(hit2, scene, recursion_depth+1);
    }
    return color;
}; //page 15


void RayTracer::RaytraceGlobal(Camera* cam, RTScene* scene, Image &image, int N) {
    std::cout << "RaytraceGlobal() called" << std::endl;
    int w = image.width; int h = image.height;
    Ray ray;    // Ray_ij,k
    glm::vec4 colorSum;     // sum of colors by N random rays at a pixel
    // c++ 11 trick to generate uniform random
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1
    float p, q;
    #pragma omp parallel for
    for (int j=0; j<h; j++){
        for (int i=0; i<w; i++){
            colorSum = glm::vec4(0.f);
            for (int k=0; k<N; k++) {
                p = dis(gen); q = dis(gen);
                ray = RandomRayThruPixel( cam, i, j, w, h, p, q);
                // pixels is 1d vector
                colorSum += FindColorRR(ray, scene);
            }
            image.pixels[(h-1-j)*w + i] = colorSum / float(N) ;
        }
    }
};

/* p,q in [0,1], generated in RaytraceGlobal() 
    Simply change those two 0.5f into p and q
*/
Ray RayTracer::RandomRayThruPixel(Camera* cam, int i, int j, int width, int height, float p, float q){
    float alpha = (float(i) + p) * 2.f/width - 1.f;
    float beta = 1.f - (float(j) + q) * 2.f/height;
    float fovy_h = 0.5f * cam->fovy / 180.f * M_PI;     // half fovy, turn degree into rad

    glm::vec3 w = glm::normalize(cam->eye - cam->target);
    glm::vec3 u = glm::normalize(glm::cross(cam->up, w));
    glm::vec3 v = glm::cross(w, u);
    glm::vec3 dir = glm::normalize(
        alpha * cam->aspect * tanf(fovy_h) * u 
        + beta * tanf(fovy_h) * v
        - w
    ); 
    return Ray(cam->eye, dir);   // return Ray in World coord
};

/* Helper function
given 2 random numbers in [0,1],
return a direction sampled from Lambert cosine hemisphere
Note this direction is under hemisphere coord (normal is [0,1,0])
 */
glm::vec3 sampleD(float s, float t) {
    float u = 2.f * M_PI * s;
    float v = sqrt(1-t*t);
    return glm::vec3(v * cos(u), t, v * sin(u));
}

/*
Goal: given d under normal [0,1,0],
    we need d' under normal hit.N
Approach: Find transformation matrix T s.t. T*[0,1,0] = hit.N
    with Rodrigues formula
Key: rotation axis is the "average" of [0,1,0] and hit.N, then rotation angle is pi
Return: d' with rotated frame
*/
glm::vec3 rotateFrame(glm::vec3 normal, glm::vec3 d) {
    // cos(pi) = -1, sin(pi) = 0
    glm::vec3 axis = glm::normalize(normal + glm::vec3(0,1,0));
    // cos(pi) * Id             (1-cos pi)
    glm::mat3 R = glm::mat3(-1.f) + 2.f * glm::outerProduct(axis, axis);
    return R * d;
}


/* Return 1/Factor x color_ij,k */
glm::vec4 RayTracer::FindColorRR(Ray ray, RTScene* scene){
    glm::vec4 color(0.f); glm::vec4 totalW(1.f);
    float factor = 1.f; int nk = 1; // path length
    Intersection hit;    // intersection with current scene;
    Material* ma;
    int terminate;  float lambda = 0.7f; // poisson probability to terminate
    glm::vec4 lightPos; glm::vec3 l;  // hit pos to light unit vec

    bool spec;  // specular or diffuse, 0.5 vs 0.5
    Ray ray2; // 2nd ray, may be mirror reflection (specular) or random (diffuse)
    float s, t;
    
    // setup uniform(0,1) generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1
    while (true) {
        hit = Intersect(ray, scene);
        if (hit.dist == INFINITY) {break;} // no hit: take bg color black
        // a random number in [0,1] is less than x with probability x
        terminate = (dis(gen) <= lambda)? 1:0;
        // Only when hit can we have triangle access
        ma = hit.triangle->material;       
        if (terminate == 1){
            factor *= lambda;
            // find diffuse color only, at hit
            glm::vec4 fiDiffuseColor(0.f);
            for (int j=0; j<scene->shader->nlights; j++) {
                lightPos = scene->shader->lightpositions[j];
                l = glm::normalize(glm::vec3(lightPos) * 1.0f - hit.P * lightPos.w);
                // Add diffuse part
                fiDiffuseColor += ma->diffuse * std::max(glm::dot(hit.N, l), 0.f) * scene->shader->lightcolors[j];
            }
            color += totalW * (ma->emision + fiDiffuseColor);            
            break;
        } else {            
            nk++;
            factor *= (1-lambda);
            color += totalW * ma->emision;            
            // (1) decide diffuse or specular, 0.5 vs 0.5
            spec = (dis(gen)<=0.5)? true:false;           
            if (spec) {   // recursive specular
                ray2 = Ray(hit.P+1e-3f*hit.N, 2.f * glm::dot(hit.N, hit.V)*hit.N - hit.V); 
                totalW *= ma->specular;                  
            } else{
                // sample a random-bounce diffuse ray; call helper functions
                s = dis(gen); t = dis(gen);
                ray2 = Ray(hit.P+1e-3f*hit.N, rotateFrame(hit.N, sampleD(s, t)) );
                totalW *= ma->diffuse;                
            }
            // recursive call, find color of next bounce
            ray = ray2;
        }
    }

    return 1/factor * color;
};

