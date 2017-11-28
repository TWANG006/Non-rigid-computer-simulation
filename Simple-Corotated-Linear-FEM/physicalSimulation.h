#ifndef _PHYSICALSIMULATION_H_
#define _PHYSICALSIMULATION_H_

#include <glm/glm.hpp>
#include "TetMesh.h"
#include "utilities.h"

#define MINIMUM(x,y) (x<y?x:y)

typedef std::map<int, glm::mat3> matrix_map;
typedef matrix_map::iterator matrix_iterator;

struct PhysicalGlobalVariables
{
	size_t iNumofTetrohedron;
	size_t total_points;
	std::vector<glm::vec3> Xi;	//Material coordinates
	std::vector<glm::vec3> P;	//world coordinates
	std::vector<glm::vec3> V;	//velocity
	std::vector<float>	MASS;	//mass matrix
	std::vector<glm::vec3> F;	//Forces
	std::vector<bool> IsFixed;		//fixed points
	glm::mat3 eye;
	
	std::vector<matrix_map> K_row;
	std::vector<matrix_map> A_row;
	std::vector<glm::vec3> F0;
	std::vector<glm::vec3> b;

	bool bUsingStiffnessWarping;

	PhysicalGlobalVariables()
		:iNumofTetrohedron(0),total_points(0),eye(glm::mat3(1)),bUsingStiffnessWarping(true)
	{}
};

extern PhysicalGlobalVariables *physicalglobals;

/* --------- Physical Quantities -----------*/
const glm::vec3 gravAcceleration = glm::vec3(0.0f, -9.81f, 0.0f);	//Gravitational acceleration	
const float fDensity = 10.0f;										//mass density
const float fDamping = 1.0f;
const float fPoisonRatio = 0.33f;			//Poisson ratio
const float fYoungModulus = 5000.0f;		//Young's modulus
const float creep = 0.20f;
const float yield = 0.04f;
const float c_max = 0.2f;

/* --------- matrix parameters used in Hookean material constitutive model -----------*/
const float fGeneralMultiplier = fYoungModulus / (1.0f + fPoisonRatio) / (1.0f - 2 * fPoisonRatio);	//d15
const float D0 = (1.0f - fPoisonRatio) * fGeneralMultiplier;									//d16
const float D1 = fPoisonRatio * fGeneralMultiplier;												//d17
const float D2 = fYoungModulus / 2 / (1.0f + fPoisonRatio);								//d18
const glm::vec3 D(D0, D1, D2);										//Matrix D isotropic elasticity

/* --------- functions for physical computing -----------*/
//void initializePhysics();	//Allocate memory for Physics
float GetTetrahedronVolume(glm::vec3 e1, glm::vec3 e2, glm::vec3 e3);	//Tetrahedron's volum = 1/6 * corresponding parrallelepiped
void addTetraheron(int i0, int i1, int i2, int i3);		//used to build the tet mesh
void genMesh(size_t xdim, size_t ydim, size_t zdim, float fWidth, float fHeight, float fDepth);		//generate the initial mesh
void calculateStiffnessK();		//Calculate stiffness matrix
void computeforce();	//Add gravity forces to nodes.
void recalcmassmatrix();	//Calculated lumped mass matrix
void stiffnessAssemble();	//Assemble the stiffness matrix
void clearStiffnessMatrix();	//Clear the stiffness matrix
glm::mat3 Gram_Schmidt(glm::mat3 G);	//Gram_Shcmidt orthoganolization
void updateOrientation();		//Compute the orientation used for warping
void resetOrientation();		//Reset orientation to I
void initPlasticity();			//Initialize the plasticity to 0
void forcePlasticity(float deltaT);				//Add plasticity
void dynamicsAssembly(float deltaT);		//Time integration
void updatePosition(float deltaT);
void collisionGround();		//Collision detection with ground
void stepPhysics(float deltaT);		//Time steps

#endif // !_PHYSICALSIMULATION_H_
