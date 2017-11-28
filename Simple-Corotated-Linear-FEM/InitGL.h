#ifndef _INITGL_H_
#define _INITGL_H_

#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/freeglut.h>
#include <glm/glm.hpp>
#include "utilities.h"

const int iGridSize = 10;

class GraphicGlobalVariables
{	
public:
	GraphicGlobalVariables()
		: iOldX(0),iOldY(0),rX(30.0f),rY(-30.0f),iState(1),fDist(-2.5f),Up(glm::vec3(0.0f,1.0f,0.0f))
	{}

	int iOldX, iOldY;
	float rX, rY;
	int iState;
	float fDist;
	glm::vec3 Up,Right,view;
	GLint iViewport[4];
	GLdouble dModelView[16];
	GLdouble dPorjection[16];

};

extern GraphicGlobalVariables *graphicalglobals;
/* --------- OpenGL initialization functions -----------*/

void initGL();
void initGlut(int argc, char **argv, char *pt_cWindowTitle, 
			  int iWindowWidth, int iWindowHeight, int &pt_iWindowId);

/* --------- Draw the environment -----------*/
void drawGround();
void drawAxis();

#endif // !_INITGL_H_


