#include "InitGL.h"

using namespace std;

GraphicGlobalVariables* graphicalglobals;

//extern void displayFunc(void);
//extern void reshapeFunc(int iWindowWidth, int iWindowHeight);
//extern void idleFunc(void);
//extern void keyboardFunc(unsigned char key, int x, int y);
//extern void onMouseDown(int button, int state, int x, int y);
//extern void onMouseMotion(int x, int y);
//extern void close();

/* --------- OpenGL initialization functions -----------*/
void initGL()
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);/*(198.0f / 256, 226.0f/256, 255.0f/256, 0.0);*/

	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glPointSize(5);
		//reshapeFunc(iWindowWidth, iWindowHeight);
	cout<<"Graphics intitialization has completed."<<endl;
}
//void initGlut(int argc, char **argv, char *pt_cWindowTitle, 
//			  int iWindowWidth, int iWindowHeight, int &pt_iWindowId)
//{
//	glutInit(&argc, argv);
//	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
//	glutInitWindowSize(iWindowWidth, iWindowHeight);
//	pt_iWindowId = glutCreateWindow(pt_cWindowTitle);
//	
//	glutDisplayFunc(displayFunc);
//	glutReshapeFunc(reshapeFunc);
//	glutIdleFunc(idleFunc);
//	glutKeyboardFunc(keyboardFunc);
//	glutMouseFunc(onMouseDown);
//	glutMotionFunc(onMouseMotion);
//	glutCloseFunc(close);
//}

/* --------- Draw the environment -----------*/
void drawGround()
{
	glBegin(GL_LINES);
		glColor3f(0.5f, 0.5f, 0.5f);
		for(int i = -iGridSize; i <= iGridSize; i++){
			glVertex3f((float)i, 0, (float)-iGridSize);
			glVertex3f((float)i, 0, (float)iGridSize);

			glVertex3f((float)-iGridSize, 0, (float)i);
			glVertex3f((float)iGridSize, 0, (float)i);
		}
	glEnd();
}
void drawAxis()
{
	glBegin(GL_LINES);
		for( int i=0; i<3; i++){
			float fColor[3] = {0.0f, 0.0f, 0.0f};
			fColor[i] = 1.0;
			glColor3fv(fColor);

			float iVertex[3] = {0.0f, 0.0f, 0.0f};
			iVertex[i] = 5.0f;
			glVertex3fv(iVertex);
			glVertex3f(0.0, 0.0, 0.0);
		}
	glEnd();
}