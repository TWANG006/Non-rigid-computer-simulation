#include "InitGL.h"
#include "physicalSimulation.h"
#include "TetMesh.h"
#include "conjugate_gradient.h"
#include "sysinfo.h"
#include <Windows.h>

#define EPSILON 0.001f
#define EPS2  EPSILON*EPSILON

using namespace std;  

const int iWidth = 1024, iHeight = 1024;

float timeStep =  1/60.0f;

static float accumulator = timeStep;
static int selected_index = -1;
vector<Tetrahedron> tetrahedra;

/* --------- functions for GLUT -----------*/
void displayFunc(void);
void reshapeFunc(int iWindowWidth, int iWindowHeight);
void idleFunc(void);
void keyboardFunc(unsigned char key, int x, int y);
void onMouseDown(int button, int state, int x, int y);
void onMouseMotion(int x, int y);
void close();
/* --------- functions for physics-----------*/
void initilization();



int main(int argc, char** argv)
{
	//initializePhysics();
	physicalglobals = new PhysicalGlobalVariables();
	//initialize_CGsolver();
	cgglobals = new ConjugateGradient();
	graphicalglobals = new GraphicGlobalVariables();
	sysinfoglobals = new Sysinfo();
	

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(iWidth, iHeight);
	glutCreateWindow("GLUT OpenCloth - Co-rotated Linear FEM with Stiffness Warping");

	glutDisplayFunc(displayFunc);
	glutReshapeFunc(reshapeFunc);
	glutIdleFunc(idleFunc);
	glutMouseFunc(onMouseDown);
	glutMotionFunc(onMouseMotion);
	glutKeyboardFunc(keyboardFunc);
	glutCloseFunc(close);

	//glewExperimental = true;
	glewInit();
	initGL();
	initilization();

	puts("Press ' ' to toggle stiffness warping on/off\n");

	glutMainLoop();

}

void onMouseDown(int button, int state, int x, int y)
{
	if( GLUT_DOWN == state){
		graphicalglobals->iOldX = x;
		graphicalglobals->iOldY = y;
		int iWindowY = (iHeight - y -1);
		int iWindowX = x;
		
		float iWindowZ = 0.0;
		glReadPixels(x, iHeight - y -1, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &iWindowZ);
		double dOldX = 0.0, dOldY =0.0, dOldZ = 0.0;
		gluUnProject(iWindowX, iWindowY,iWindowZ, 
			graphicalglobals->dModelView,graphicalglobals->dPorjection,
			graphicalglobals->iViewport,&dOldX, &dOldY,&dOldZ);
		glm::vec3 pt(dOldX, dOldY, dOldZ);

		//printf_s("\nObject [ %3.3f, %3,3f, %3,3f ]",dOldX,dOldY,dOldZ);
		for(size_t i=0; i<physicalglobals->total_points; i++){
			if(glm::distance(physicalglobals->P[i],pt)<0.01){
				selected_index =i;

				printf_s("Point %d is picked:\n",i);
				printf_s("Pt [%3.3f, %3.3f, %3.3f]\n",physicalglobals->P[i].x,physicalglobals->P[i].y,physicalglobals->P[i].z);
				break;	
			}
		}
	}
	if(GLUT_MIDDLE_BUTTON == button)
		graphicalglobals->iState = 0;
	else
		graphicalglobals->iState = 1;
	if(GLUT_UP == state){
		selected_index = -1;
		updateOrientation();
		glutSetCursor(GLUT_CURSOR_INHERIT);
	}
}
void onMouseMotion(int x, int y)
{
	if(selected_index == -1){
		if(graphicalglobals->iState == 0)
			graphicalglobals->fDist *= (1.0 + (y-graphicalglobals->iOldY)/60.0);
		else{
			graphicalglobals->rY +=(x - graphicalglobals->iOldX)/5.0f;
			graphicalglobals->rX +=(y - graphicalglobals->iOldY)/5.0f;
		}
	}
	else{
		float fDelta = 1000/abs(graphicalglobals->fDist);
		float fValX = (x - graphicalglobals->iOldX)/fDelta;
		float fValY = (graphicalglobals->iOldY - y)/fDelta;
		if(abs(fValX)>abs(fValY))
			glutSetCursor(GLUT_CURSOR_LEFT_RIGHT);
		else
			glutSetCursor(GLUT_CURSOR_UP_DOWN);
		physicalglobals->V[selected_index] = glm::vec3(0);
		physicalglobals->P[selected_index].x += graphicalglobals->Right[0]*fValX;
		float fNewval = physicalglobals->P[selected_index].y + graphicalglobals->Up[1]*fValY;
		if(fNewval >0)
			physicalglobals->P[selected_index].y = fNewval;
		physicalglobals->P[selected_index].z += graphicalglobals->Right[2]*fValX + graphicalglobals->Up[2] * fValY;
	}
	graphicalglobals->iOldX = x;
	graphicalglobals->iOldY = y;

	glutPostRedisplay();
}
void keyboardFunc(unsigned char key, int x, int y)
{
	switch(key){
	case ' ': 
		physicalglobals->bUsingStiffnessWarping = !physicalglobals->bUsingStiffnessWarping;
		break;
	case 27:
		exit(0);
	}
	printf("Stiffness Warping %s\n", physicalglobals->bUsingStiffnessWarping?"On":"Off");
	glutPostRedisplay();
}
void initilization()
{
	genMesh(10,3,3,0.1f,0.1f,0.1f);
	physicalglobals->iNumofTetrohedron = tetrahedra.size();

	cout<<physicalglobals->iNumofTetrohedron<<tetrahedra.size();

	physicalglobals->total_points = physicalglobals->P.size();
	physicalglobals->MASS.resize(physicalglobals->total_points);

	physicalglobals->A_row.resize(physicalglobals->total_points);
	physicalglobals->K_row.resize(physicalglobals->total_points);
	physicalglobals->b.resize(physicalglobals->total_points);
	physicalglobals->V.resize(physicalglobals->total_points);
	physicalglobals->F.resize(physicalglobals->total_points);
	physicalglobals->F0.resize(physicalglobals->total_points);
	cgglobals->residual.resize(physicalglobals->total_points);
	cgglobals->q.resize(physicalglobals->total_points);
	cgglobals->d.resize(physicalglobals->total_points);

	memset(&(physicalglobals->V[0].x),0,physicalglobals->total_points*sizeof(glm::vec3));

	sysinfoglobals->fStartTime = (float)glutGet(GLUT_ELAPSED_TIME);
	sysinfoglobals->currentTime = sysinfoglobals->fStartTime;

	QueryPerformanceFrequency(&sysinfoglobals->frequency);

	QueryPerformanceCounter(&sysinfoglobals->t1);

	wglSwapIntervalEXT(0);

	calculateStiffnessK();
	clearStiffnessMatrix();
	recalcmassmatrix();
	initPlasticity();
	
}
void idleFunc(void)
{
	if(accumulator >= timeStep){
		stepPhysics(timeStep);
		accumulator -= timeStep;
	}
	glutPostRedisplay();
}
void reshapeFunc(int iWindowWidth, int iWindowHeight)
{
	glViewport(0, 0, iWindowWidth, iWindowHeight);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, (GLfloat)iWindowWidth / (GLfloat)iWindowHeight, 0.1f, 100.0f);

	glGetIntegerv(GL_VIEWPORT, graphicalglobals->iViewport);
	glGetDoublev(GL_PROJECTION_MATRIX, graphicalglobals->dPorjection);

	glMatrixMode(GL_MODELVIEW);
}
void displayFunc(void)
{
	size_t i=0;
	float newtime = (float)glutGet(GLUT_ELAPSED_TIME);
	sysinfoglobals->fFrameTime = newtime - sysinfoglobals->currentTime;
	sysinfoglobals->currentTime = newtime;

	QueryPerformanceCounter(&sysinfoglobals->t2);
	sysinfoglobals->fFrameTimeQP = (sysinfoglobals->t2.QuadPart - sysinfoglobals->t1.QuadPart) * 1000.0f / sysinfoglobals->frequency.QuadPart;
	sysinfoglobals->t1 = sysinfoglobals->t2;
	accumulator += sysinfoglobals->fFrameTimeQP;

	++sysinfoglobals->iTotalFrames;
	if((newtime - sysinfoglobals->fStartTime) > 1000.00)
	{
		float elapsedTime = (newtime - sysinfoglobals->fStartTime);
		sysinfoglobals->fps = (sysinfoglobals->iTotalFrames / elapsedTime) * 1000.0;
		sysinfoglobals->fStartTime = newtime;
		sysinfoglobals->iTotalFrames = 0;
	}

	sprintf_s(sysinfoglobals->cTitleinfo, "FPS: %3.2f, Frame time (GLUT): %3.4f msecs, Frame time (QP): %3.3f, Stiffness Warping: %s", sysinfoglobals->fps, sysinfoglobals->fFrameTime, sysinfoglobals->fFrameTimeQP,physicalglobals->bUsingStiffnessWarping?"On":"Off");
	glutSetWindowTitle(sysinfoglobals->cTitleinfo);
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0.96f, 0.96f, 0.96f, 0.0f);
	glLoadIdentity();
	glTranslatef(0.0, 0.0, graphicalglobals->fDist);
	glRotatef(graphicalglobals->rX, 1.0, 0.0, 0.0);
	glRotatef(graphicalglobals->rY, 0.0, 1.0, 0.0);

	glGetDoublev(GL_MODELVIEW_MATRIX, graphicalglobals->dModelView);
	graphicalglobals->view.x = (float)-graphicalglobals->dModelView[2];
	graphicalglobals->view.y = (float)-graphicalglobals->dModelView[6];
	graphicalglobals->view.z = (float)-graphicalglobals->dModelView[10];
	graphicalglobals->Right = glm::cross(graphicalglobals->view,graphicalglobals->Up);

	//drawGround();
	glLineWidth(4.0);
	drawAxis();

	glLineWidth(1.0);
	glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_LINES);
		for(i = 0; i <physicalglobals->iNumofTetrohedron;i++){
			int i0 = tetrahedra[i].iIndex[0];
			int i1 = tetrahedra[i].iIndex[1];
			int i2 = tetrahedra[i].iIndex[2];
			int i3 = tetrahedra[i].iIndex[3];
			glm::vec3 p0 = physicalglobals->P[i0];
			glm::vec3 p1 = physicalglobals->P[i1];
			glm::vec3 p2 = physicalglobals->P[i2];
			glm::vec3 p3 = physicalglobals->P[i3];

			glVertex3f(p3.x, p3.y, p3.z);		glVertex3f(p0.x, p0.y, p0.z);
			glVertex3f(p3.x, p3.y, p3.z);		glVertex3f(p1.x, p1.y, p1.z);
			glVertex3f(p3.x, p3.y, p3.z);		glVertex3f(p2.x, p2.y, p2.z);

			glVertex3f(p0.x, p0.y, p0.z);		glVertex3f(p1.x, p1.y, p1.z);
			glVertex3f(p0.x, p0.y, p0.z);		glVertex3f(p2.x, p2.y, p2.z);

			glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p2.x, p2.y, p2.z);
		}
	glEnd();

	glBegin(GL_POINTS);
		for(i = 0; i< physicalglobals->total_points; i++){
			glm::vec3 p = physicalglobals->P[i];
			int is = (i == selected_index);
			glColor3f((float)!is, (float)is, (float)is);
			glVertex3f(p.x, p.y, p.z);
		}
	glEnd();
	
	glutSwapBuffers();
}
void close()
{
	delete physicalglobals;
	delete graphicalglobals;
	delete cgglobals;
	delete sysinfoglobals;
}