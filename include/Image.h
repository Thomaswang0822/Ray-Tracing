
#ifndef __IMAGE_H__
#define __IMAGE_H__
#include <vector>
#include <iostream>
#include <glm/glm.hpp>
#include <OpenGL/gl3.h>
#include <OpenGL/glext.h>
#include <GLUT/glut.h>

class Image {
private:
    unsigned int fbo; // frame buffer object
    unsigned int tbo; // texture buffer object

public: 
    int width;
    int height;
    std::vector<glm::vec3> pixels; // RGB colors

    void init(int w, int h) {
        width = w;
        height = h;
        glGenFramebuffers(1, &fbo);
        glGenTextures(1, &tbo);
    }

    void draw(void) {
        glBindTexture(GL_TEXTURE_2D, tbo);
        // load texture pointing to start of RGB vector
        glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,width,height,
                    0, GL_RGB, GL_FLOAT, &pixels[0][0]); 
        glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
        // attach texture and the read frame
        glFramebufferTexture2D(GL_READ_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                    GL_TEXTURE_2D, tbo, 0);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0); // if not already so
        // copy framebuffer from read to write
        glBlitFramebuffer(0,0,width,height,0,0,width,height,
                    GL_COLOR_BUFFER_BIT, GL_NEAREST);
    }

    // Helper function: populate pixels
    void fillPixels() {
        int length = width * height;
        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                // Give the dummy photo transition of colors
                glm::vec3 pix = glm::vec3(i/width, j/height, i*j/length);
                pixels.push_back(pix);
            }
        }
    }
};


#endif