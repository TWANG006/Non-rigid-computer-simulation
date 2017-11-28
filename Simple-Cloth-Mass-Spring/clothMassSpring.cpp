#include <gl/glew.h>
#include <gl/wglew.h>
#include <gl/freeglut.h>
#include <vector>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "Spring.h"

using namespace std;

const int WIDTH = 1024, HEIGHT= 1024;

//Initial values needed to initilize the regular grid M*N
int inumM = 20, inumN = 20;
const size_t stTotalpoints = (inumM+1) * (inumN+1);
int isize =4;
float fhsize = 4/2.0f;

//Set timesteps for the animation
float fTimestep = 1/60.0f;
float fCurrentTime = 0.0;
double dAccumulator = fTimestep;
int iSelectedPt = -1;

//Initial values for springs
vector<GLushort> indices;
vector<Spring> springs;
vector<glm::vec3> X;
vector<glm::vec3> V;
vector<glm::vec3> F;
vector<glm::vec3> sumF;  
vector<glm::vec3> sumV;

//initial values for the initial state
int iOldX = 0, iOldY = 0;
float fRotationX = 15, fRotationY = 0;
int iState = 1;
float fDist = -10;
const int GROUND_SIZE = 10;

const int STRUCTURAL_SPRING = 0;
const int SHEAR_SPRING = 1;
const int BEND_SPRING = 2;
int iSpringcount = 0;

//Initial Spring parameters
const float DEFAULT_DAMPING = -0.0125f;
float fKsStretch = 0.5f, fKdStretch = -0.25f;
float fKsShear = 0.5f, fKdShear = -0.25f;
float fKsBending = 0.85f, fKdBending = -0.25f;
glm::vec3 gravity(0.0f, -0.00981f, 0.0f);
//glm::vec3 wind(0.005f,0.0f,0.0f);
const float MASS = 0.5f;

//GL initialization
GLint viewport[4];
GLdouble ModelView[16];
GLdouble Projection[16];
glm::vec3 Up(0,1,0);
glm::vec3 viwDir, Right;


LARGE_INTEGER frequency;
LARGE_INTEGER t1, t2;
double dFrameTimeQP = 0;
float fFrameTime = 0;

//initialize the displayed events
float fStarttime = 0, fps = 0;
int iTotalframes = 0;

char info[MAX_PATH] = {0};

//Set the barriers to ellipsoid
glm::mat4 ellipsoid, inverse_ellipsoid;
int iStacks = 30;
int iSlices = 30;
float fRadius = 1;
glm::vec3 center(0,0,0);
float radius = 1;

void StepPhysics(float dt);

//Begin to implement
void AddingSpring(int iP1, int iP2, float fKs, float fKd, int iType)
{
	Spring spring;
	spring.iP1 = iP1;
	spring.iP2 = iP2;
	spring.fKs = fKs;
	spring.fKd = fKd;
	spring.iType = iType;
	glm::vec3 deltaP = X[iP1] - X[iP2];
	spring.fRestlength = sqrt(glm::dot(deltaP,deltaP));
	springs.push_back(spring);
}

//Mouse and key Events
void OnMouse(int button, int state, int x, int y)
{
	if(GLUT_DOWN == state)
	{
		iOldX = x, iOldY =y;
		int iWindowX = x;
		int iWindowY = HEIGHT - y;
		float fNormalizedX = float(x)/float(WIDTH/2.0);
		float fNormalizedY = float(HEIGHT-y)/float(HEIGHT/2.0);
		float iZ = 0;
		glReadPixels(x,HEIGHT-y,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,&iZ);
		if(1 == iZ)
			iZ = 0;
		double dObjectX = 0, dObjectY =0, dObjectZ = 0;
		gluUnProject(iWindowX, iWindowY,iZ, ModelView,Projection,viewport,&dObjectX, &dObjectY, &dObjectZ);
		glm::vec3 pt(dObjectX, dObjectY, dObjectZ);
		for(size_t i = 0; i<stTotalpoints; i++)
		{
			if(glm::distance(X[i],pt)<0.1)
			{
				iSelectedPt = i;
				cout<<"Point picked at position"<<i<<endl;
				break;
			}
		}
	}
	if(GLUT_MIDDLE_BUTTON == button)
		iState = 0;
	else
		iState = 1;
	if(GLUT_UP == state)
	{
		iSelectedPt = -1;
		glutSetCursor(GLUT_CURSOR_INHERIT);
	}
}
void Motion(int x, int y)
{
	if(-1 == iSelectedPt)
	{
		if(iState == 0)
			fDist *= (1 + float(y - iOldY)/60.0f);
		else
		{
			fRotationY += float(x - iOldX)/5.0f;
			fRotationX += float(y - iOldY)/5.0f;
		}
	}
	else
	{
		float fDelta = 1500/abs(fDist);
		float fValX = (x - iOldX)/fDelta;
		float fValY = (iOldY - y)/fDelta;
		if(abs(fValX) > abs(fValY))
			glutSetCursor(GLUT_CURSOR_LEFT_RIGHT);
		else
			glutSetCursor(GLUT_CURSOR_UP_DOWN);

		V[iSelectedPt] = glm::vec3(0);
		X[iSelectedPt].x += Right[0]*fValX;
		float fNewvalue = X[iSelectedPt].y+Up[1]*fValY;
		if(fNewvalue >0)
			X[iSelectedPt].y = fNewvalue;
		X[iSelectedPt].z += Right[2]*fValX + Up[2]*fValY;

		V[iSelectedPt].x = 0;
		V[iSelectedPt].y = 0;
		V[iSelectedPt].z = 0;
	}
	iOldX =x;
	iOldY =y;

	glutPostRedisplay();
}
void keyPress(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 27:
		  {exit(0);}
		  break;
	default:
		break;
	}
}


//Function to initialize OpenGL
void DrawGround()
{
	glBegin(GL_LINES);
		glColor3f(1.0,1.0,1.0);
		for(int i = -GROUND_SIZE; i<=GROUND_SIZE; i++)
		{
			glVertex3f((float)i,0,(float)-GROUND_SIZE);
			glVertex3f((float)i,0,(float)GROUND_SIZE);

			glVertex3f((float)-GROUND_SIZE,0,(float)i);
			glVertex3f((float)GROUND_SIZE,0,(float)i);
		}
	glEnd();
}
void Init()
{
	fStarttime = (float)glutGet(GLUT_ELAPSED_TIME);
	fCurrentTime =  fStarttime;

	//ticks per second
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&t1);

	glEnable(GL_DEPTH_TEST);
	int count=0;
	int l1=0, l2=0;
	float fYpos = 7.0f;
	int m = inumM+1, n = inumN+1;

	indices.resize(inumM*inumN*2*3);
	X.resize(stTotalpoints);
	V.resize(stTotalpoints);
	F.resize(stTotalpoints);

	sumF.resize(stTotalpoints);
	sumV.resize(stTotalpoints);

	//Bulid up positions X
	for(int j=0; j<n; j++)
	{
		for(int i=0; i<m; i++)
		{
			X[count++] = glm::vec3( ((float(i)/(m-1))*2-1)*fhsize,isize+1,((float(j)/(n-1))*isize));
			//X.push_back(glm::vec3( ((float(i)/(m-1))*2-1)*fhsize,isize+1,((float(j)/(n-1))*isize)));
		}
	}
	//Build up velocities V
	memset(&(V[0].x),0,stTotalpoints*sizeof(glm::vec3));
	//Bulid up Indices
	GLushort* id = &indices[0];
	for (int i=0; i<inumN; i++)
	{
		for(int j=0; j<inumM; j++)
		{
			int i0 = i*(inumM+1) + j;
			int i1 = i0 + 1;
			int i2 = i0 + (inumM+1);
			int i3 = i2+1;
			if((j+i)%2)
			{
				*id++ = i0; *id++ = i2; *id++ = i1;
				*id++ = i1; *id++ = i2; *id++ = i3;
			}
			else
			{
				*id++ = i0; *id++ = i2; *id++ = i3;
				*id++ = i0; *id++ = i3; *id++ = i1;
			}
		}
	}
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);   // Can be changed
	glPointSize(5);

	wglSwapIntervalEXT(0); //disable vsync

	//Build up springs
	//horizontal and vertical
	for(l1 = 0; l1<n; l1++)
		for(l2=0; l2<(m-1); l2++)
		{
			AddingSpring((l1*m)+l2,(l1*m)+l2+1,fKsStretch,fKdStretch,STRUCTURAL_SPRING);
		}
	for(l1 = 0; l1<m ;l1++)
		for(l2 = 0; l2<(n-1); l2++)
		{
			AddingSpring((l2*m)+l1,((l2+1)*m) +l1,fKsStretch,fKdStretch,STRUCTURAL_SPRING);
		}
	//Shearing springs
	for(l1 = 0; l1< n-1; l1++)
	{
		for(l2=0; l2< m-1; l2++)
		{
			AddingSpring((l1*m)+l2,((l1+1)*m)+l2+1, fKsShear,fKdShear, SHEAR_SPRING);
			AddingSpring(((l1+1)*m)+l2, (l1*m)+l2+1,fKsShear,fKdShear, SHEAR_SPRING);
		}
	}
	//Bending springs
	for(l1=0; l1<n; l1++)
	{
		for(l2=0; l2<m-2; l2++)
		{
			AddingSpring((l1*m)+l2,(l1*m)+l2+2,fKsBending,fKdBending,BEND_SPRING);
		}
		AddingSpring((l1*m)+(m-3), (l1*m)+(m-1),fKsBending,fKdBending,BEND_SPRING);
	}
		
	for(l1=0; l1<m; l1++)
	{
		for( l2=0; l2<n-2; l2++)
		{
			AddingSpring((l2*m)+l1, ((l2+2)*m)+l1,fKsBending,fKdBending,BEND_SPRING);	
		}
		AddingSpring(((n-3)*m)+l1,((n-1)*m)+l1,fKsBending, fKdBending, BEND_SPRING);
	}
}
void reshape(int w, int h)
{
	glViewport(0,0,w,h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, (GLfloat)w/(GLfloat)h,1.0f, 100.0f);

	glGetIntegerv(GL_VIEWPORT,viewport);
	glGetDoublev(GL_PROJECTION_MATRIX,Projection);

	glMatrixMode(GL_MODELVIEW);
}
void display()
{
	size_t i = 0;
	float fNewTime = (float)glutGet(GLUT_ELAPSED_TIME);
	fFrameTime = fNewTime - fCurrentTime;
	fCurrentTime = fNewTime;

	QueryPerformanceCounter(&t2);

	dFrameTimeQP = (t2.QuadPart - t1.QuadPart) * 1000.0/frequency.QuadPart;
	t1=t2;
	dAccumulator += dFrameTimeQP;

	++iTotalframes;
	if((fNewTime-fStarttime)>1000)
	{
		float fElapsedTime = fNewTime-fStarttime;
		fps = (iTotalframes / fElapsedTime) * 1000;
		fStarttime = fNewTime;
		iTotalframes = 0;
	}

	sprintf_s(info, "FPS: %3.2f, Frame time (GLUT): %3.5f msecs, Frame time (QP): %3.3f", fps, fFrameTime,dFrameTimeQP);
	glutSetWindowTitle(info);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0.5,0.5,0.5,1.0);
	glLoadIdentity();
	glTranslatef(0,0,fDist);
	glRotatef(fRotationX,1,0,0);
	glRotatef(fRotationY,0,1,0);

	glGetDoublev(GL_MODELVIEW_MATRIX,ModelView);
	viwDir.x = (float)-ModelView[2];
	viwDir.y = (float)-ModelView[6];
	viwDir.z = (float)-ModelView[10];
	Right = glm::cross(viwDir,Up);

	DrawGround();

	//Draw ellipsoid
	/*glColor3f(0,1,0);
	glPushMatrix();
		glMultMatrixf(glm::value_ptr(ellipsoid));
			glutWireSphere(fRadius, iSlices, iStacks);
	glPopMatrix();*/

	//Draw poylgons
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glColor3f(0,0.5,0.25);
	glBegin(GL_TRIANGLES);
		for(i = 0; i<indices.size(); i+=3)
		{
			glm::vec3 p1 = X[indices[i]];
			glm::vec3 p2 = X[indices[i+1]];
			glm::vec3 p3 = X[indices[i+2]];
			glVertex3f(p1.x,p1.y,p1.z);
			glVertex3f(p2.x,p2.y,p2.z);
			glVertex3f(p3.x,p3.y,p3.z);
		}
	glEnd();

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor3f(1,0.5,0.25);
	glBegin(GL_TRIANGLES);
		for(i = 0; i<indices.size(); i+=3)
		{
			glm::vec3 p1 = X[indices[i]];
			glm::vec3 p2 = X[indices[i+1]];
			glm::vec3 p3 = X[indices[i+2]];
			glVertex3f(p1.x,p1.y,p1.z);
			glVertex3f(p2.x,p2.y,p2.z);
			glVertex3f(p3.x,p3.y,p3.z);
		}
	glEnd();

	//Draw particles
	glBegin(GL_POINTS);
		for(i = 0; i<stTotalpoints;i++)
		{
			glm::vec3 p = X[i];
			int is = (i==iSelectedPt);
			glColor3f((float)!is,(float)!is,(float)is);
			glVertex3f(p.x,p.y,p.z);
		}
	glEnd();

	glutSwapBuffers();
}
void shutdown()
{
	X.clear();
	V.clear();
	F.clear();
	sumF.clear();
	sumV.clear();
	indices.clear();
	springs.clear();
}


//
//Physical simulation and Implementation
//Compute External forces
void ComputeForces()
{
	size_t i=0;
	//Add gravity and wind
	for(i = 0; i<stTotalpoints;i++)
	{
		F[i] = glm::vec3(0);

		//fix both ends of the first row
		if(i!=0 && i!=inumM)
		{
			F[i] += gravity;
			//F[i] += wind;
		}
			
		F[i] += DEFAULT_DAMPING*V[i];
	}
	//Add spring forces
	for(i = 0; i<springs.size(); i++)
	{
		glm::vec3 p1 = X[springs[i].iP1];
		glm::vec3 p2 = X[springs[i].iP2];
		glm::vec3 v1 = V[springs[i].iP1];
		glm::vec3 v2 = V[springs[i].iP2];
		glm::vec3 deltaP = p1 - p2;
		glm::vec3 deltaV = v1 - v2;
		float dist = glm::length(deltaP);

		float f1 = -springs[i].fKs * (dist-springs[i].fRestlength);
		float f2 = springs[i].fKd * (glm::dot(deltaV,deltaP)/dist);
		glm::vec3 springForce = (f1+f2) * glm::normalize(deltaP);

		if(springs[i].iP1 != 0 && springs[i].iP1 != inumM)
			F[springs[i].iP1] += springForce;
		if(springs[i].iP2 != 0 && springs[i].iP2 != inumM)
			F[springs[i].iP2] -= springForce;
	}
}


//Explicit Euler solver
void ExplicitEuler(float deltaT)
{
	float deltaTimeMass = deltaT/MASS;
	size_t i = 0;
	for(i=0; i<stTotalpoints; i++)
	{
		glm::vec3 oldV = V[i];
		V[i] += (F[i]*deltaTimeMass);
		X[i] += deltaT*V[i];
		
		if(X[i].y<0)
			X[i].y = 0;
	}
}
void ApplyProvotDynamicInverse()
{
	for(size_t i=0; i<springs.size(); i++)
	{
		glm::vec3 p1 = X[springs[i].iP1];
		glm::vec3 p2 = X[springs[i].iP2];
		glm::vec3 deltaP = p1-p2;
		float dist = glm::length(deltaP);
		if(dist>springs[i].fRestlength)
		{
			dist -= (springs[i].fRestlength);
			dist /=2.0f;
			deltaP = glm::normalize(deltaP);
			deltaP *= dist;
			if(springs[i].iP1 == 0 || springs[i].iP1 == inumM)
			{
				V[springs[i].iP2] += deltaP;
			}
			else if(springs[i].iP2 ==0 || springs[i].iP2 == inumM)
			{
				V[springs[i].iP1] -= deltaP;
			}
			else
			{
				V[springs[i].iP1] -= deltaP;
				V[springs[i].iP2] += deltaP;
			}
		}
	}
}
void idle()
{
	if (dAccumulator >=fTimestep)
	{
		StepPhysics(fTimestep);
		dAccumulator -= fTimestep;
	}
	glutPostRedisplay();
}
void StepPhysics(float dt)
{
	ComputeForces();
	ExplicitEuler(dt);
	/*ellipsoidCollision();*/
	ApplyProvotDynamicInverse();
}

void main( int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA |GLUT_DEPTH);
	glutInitWindowSize(WIDTH,HEIGHT);
	glutCreateWindow("Explicit Euler Surface Simulation");

	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutIdleFunc(idle);

	glutMouseFunc(OnMouse);
	glutMotionFunc(Motion);
	glutKeyboardFunc(keyPress);

	glutCloseFunc(shutdown);
	

	glewInit();
	Init();

	glutMainLoop();

}