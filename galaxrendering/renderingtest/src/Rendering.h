/*
 * Rendering.h
 *
 *  Created on: May 14, 2015
 *      Author: sarwan
 */

#ifndef RENDERING_H_
#define RENDERING_H_

#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glext.h>
#include <GL/glu.h>
#include <SDL/SDL.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include<string.h>
#include<time.h>

struct prop {
	float x;
	float y;
	float z;
	float mass;
	float xcol;
	float ycol;
	float zcol;
};

void initGL();

void initGLFW();
void initPointSpriteExt();
static void error_callback(int error, const char* description);
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void Pre_Render();
void Render(std::vector<prop> coord, int num);
void Post_Render();
//void SaveToTGA(int idx=-1);  	not yet working : taking snapshots
//void SaveToTGA2(const std::string &sName);		not yet working

extern GLuint m_texStar;
extern GLFWwindow* window;
//extern int m_idxSnapshot;




#endif /* RENDERING_H_ */
