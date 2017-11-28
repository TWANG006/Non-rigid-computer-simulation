#include "physicalSimulation.h"
#include "conjugate_gradient.h"
#include <limits>

using namespace std;

//Initialize the physical parameters
/*void initializePhysics()
{
	physicalglobals = new PhysicalGlobalVariables();
}*/
//compute the volume of a tetraheron

PhysicalGlobalVariables* physicalglobals;

float GetTetrahedronVolume(glm::vec3 e1, glm::vec3 e2, glm::vec3 e3)
{
	return (glm::dot(e1, glm::cross(e2,e3)))/6.0f;
}
//Add one new tetraheron to the existing tetmesh vector
void addTetraheron(int i0, int i1, int i2, int i3)
{
	Tetrahedron t;
	t.iIndex[0] = i0;
	t.iIndex[1] = i1;
	t.iIndex[2] = i2;
	t.iIndex[3] = i3;

	tetrahedra.push_back(t);
}
//Generate the initial mesh
void genMesh(size_t xdim, size_t ydim, size_t zdim, float fWidth, float fHeight, float fDepth)
{
	physicalglobals->total_points = (xdim+1)*(ydim+1)*(zdim+1);
	physicalglobals->P.resize(physicalglobals->total_points);
	physicalglobals->Xi.resize(physicalglobals->total_points);
	physicalglobals->IsFixed.resize(physicalglobals->total_points);
	int ind = 0;
	
	//build every tetrohedron's vertex
	for(unsigned int  x=0; x<=xdim; x++){
		for(unsigned int y=0; y<=ydim; y++){
			for(unsigned int z=0; z<=zdim; z++){
				physicalglobals->Xi[ind] = physicalglobals->P[ind] 
				= glm::vec3(fWidth*x,fHeight*z,fDepth*y);
				if(physicalglobals->Xi[ind].x < 0.01)
					physicalglobals->IsFixed[ind] = true;

				ind++;
			}
		}
	}
	for(size_t i=0; i<physicalglobals->total_points; i++){
		physicalglobals->P[i].y += 0.5;
		physicalglobals->P[i].z -= zdim/2 * fDepth;
	}
	for (size_t i = 0; i < xdim; ++i) {
		for (size_t j = 0; j < ydim; ++j) {
			for (size_t k = 0; k < zdim; ++k) {
				int p0 = (i * (ydim + 1) + j) * (zdim + 1) + k;
				int p1 = p0 + 1;
				int p3 = ((i + 1) * (ydim + 1) + j) * (zdim + 1) + k;
				int p2 = p3 + 1;
				int p7 = ((i + 1) * (ydim + 1) + (j + 1)) * (zdim + 1) + k;
				int p6 = p7 + 1;
				int p4 = (i * (ydim + 1) + (j + 1)) * (zdim + 1) + k;
				int p5 = p4 + 1;
				// Ensure that neighboring tetras are sharing faces
				if ((i + j + k) % 2 == 1) {
					addTetraheron(p1,p2,p6,p3);
					addTetraheron(p3,p6,p4,p7);
					addTetraheron(p1,p4,p6,p5);
					addTetraheron(p1,p3,p4,p0);
					addTetraheron(p1,p6,p4,p3); 
				} else {
					addTetraheron(p2,p0,p5,p1);
					addTetraheron(p2,p7,p0,p3);
					addTetraheron(p2,p5,p7,p6);
					addTetraheron(p0,p7,p5,p4);
					addTetraheron(p2,p0,p7,p5); 
				}
				physicalglobals->iNumofTetrohedron += 5;
			}
		}
	}
}
//Calculate Stiffness K
void calculateStiffnessK()
{
	for(size_t k=0; k<tetrahedra.size(); k++){
		glm::vec3 x0 = physicalglobals->Xi[tetrahedra[k].iIndex[0]];
		glm::vec3 x1 = physicalglobals->Xi[tetrahedra[k].iIndex[1]];
		glm::vec3 x2 = physicalglobals->Xi[tetrahedra[k].iIndex[2]];
		glm::vec3 x3 = physicalglobals->Xi[tetrahedra[k].iIndex[3]];

		glm::vec3 e10 = x1-x0;
		glm::vec3 e20 = x2-x0;
		glm::vec3 e30 = x3-x0;

		tetrahedra[k].e1 = e10;
		tetrahedra[k].e2 = e20;
		tetrahedra[k].e3 = e30;

		tetrahedra[k].fVolume = GetTetrahedronVolume(e10, e20, e30);

		glm::mat3 E = glm::mat3(e10.x, e10.y, e10.z,
							    e20.x, e20.y, e20.z,
								e30.x, e30.y, e30.z);

		float detE = glm::determinant(E);
		float invDetE = 1.0f / detE;

		float invE10 = (e20.z*e30.y - e20.y*e30.z)*invDetE;
		float invE20 = (e10.y*e30.z - e10.z*e30.y)*invDetE;
		float invE30 = (e10.z*e20.y - e10.y*e20.z)*invDetE;
		float invE00 = -invE10-invE20-invE30;

		float invE11 = (e20.x*e30.z - e20.z*e30.x)*invDetE;
		float invE21 = (e10.z*e30.x - e10.x*e30.z)*invDetE;
		float invE31 = (e10.x*e20.z - e10.z*e20.x)*invDetE;
		float invE01 = -invE11-invE21-invE31;

		float invE12 = (e20.y*e30.x - e20.x*e30.y)*invDetE;
		float invE22 = (e10.x*e30.y - e10.y*e30.x)*invDetE;
		float invE32 = (e10.y*e20.x - e10.x*e20.y)*invDetE;
		float invE02 = -invE12-invE22-invE32;

		tetrahedra[k].B[0] = glm::vec3(invE00, invE01, invE02);
		tetrahedra[k].B[1] = glm::vec3(invE10, invE11, invE12);
		tetrahedra[k].B[2] = glm::vec3(invE20, invE21, invE22);
		tetrahedra[k].B[3] = glm::vec3(invE30, invE31, invE32);

		for(int i=0; i< 4; i++){
			for(int j=0; j< 4; j++){
				glm::mat3& Ke = tetrahedra[k].Ke[i][j];
				float bn = tetrahedra[k].B[i].x;
				float cn = tetrahedra[k].B[i].y;
				float dn = tetrahedra[k].B[i].z;
				float bm = tetrahedra[k].B[j].x;
				float cm = tetrahedra[k].B[j].y;
				float dm = tetrahedra[k].B[j].z;

				Ke[0][0] = D0*bn*bm + D2*(cn*cm+dn*dm);
				Ke[0][1] = D1*bn*cm + D2*cn*bm;
				Ke[0][2] = D1*bn*dm + D2*dn*bm;
				Ke[1][0] = D1*cn*bm + D2*bn*cm;
				Ke[1][1] = D0*cn*cm + D2*(bn*bm+dn*dm);
				Ke[1][2] = D1*cn*dm + D2*dn*cm;
				Ke[2][0] = D1*dn*bm + D2*bn*dm;
				Ke[2][1] = D1*dn*cm + D2*cn*dm;
				Ke[2][2] = D0*dn*dm + D2*(bn*bm+cn*cm);

				Ke *= tetrahedra[k].fVolume;
			}
		}
	}
}
//Compute vetex force
void computeforce()
{
	size_t i=0;
	for(i = 0; i< physicalglobals->total_points; i++){
		physicalglobals->F[i] = glm::vec3(0);

		//Only consider the gravity for non-fixed points
		physicalglobals->F[i] += gravAcceleration*physicalglobals->MASS[i];
	}
}
//Reassemble the lumped mass matrix
void recalcmassmatrix()
{
	for(size_t i=0; i< physicalglobals->total_points;i++){
		if(physicalglobals->IsFixed[i])
			physicalglobals->MASS[i] = numeric_limits<float>::max( );
		else
			physicalglobals->MASS[i] = 1.0f/physicalglobals->total_points;
	}
	for(size_t i =0; i < physicalglobals->iNumofTetrohedron; i++){
		float fM = (fDensity * tetrahedra[i].fVolume)*0.25f;
		physicalglobals->MASS[tetrahedra[i].iIndex[0]] += fM;
		physicalglobals->MASS[tetrahedra[i].iIndex[1]] += fM;
		physicalglobals->MASS[tetrahedra[i].iIndex[2]] += fM;
		physicalglobals->MASS[tetrahedra[i].iIndex[3]] += fM;
	}
}
//Assemble the stiffness matrix K
void stiffnessAssemble()
{
	for (size_t e = 0; e < physicalglobals->iNumofTetrohedron; e++){
		glm::mat3 Re = tetrahedra[e].Re;
		glm::mat3 ReT = glm::transpose(Re);

		for(int i = 0; i < 4; i++){
			glm::vec3 vTempForce = glm::vec3(0.0f,0.0f,0.0f);
			for(int j=0; j < 4; j++){
				glm::mat3 tmpKe = tetrahedra[e].Ke[i][j];
				glm::vec3 x0 = physicalglobals->Xi[tetrahedra[e].iIndex[j]];
				glm::vec3 prod = glm::vec3(tmpKe[0][0]*x0.x + tmpKe[0][1]*x0.y + tmpKe[0][2]*x0.z,
										   tmpKe[1][0]*x0.x + tmpKe[1][1]*x0.y + tmpKe[1][2]*x0.z,
										   tmpKe[2][0]*x0.x + tmpKe[2][1]*x0.y + tmpKe[2][2]*x0.z);
				vTempForce += prod;
				if ( j >= i){
					glm::mat3 tmp = Re * tmpKe * ReT;
					int index = tetrahedra[e].iIndex[i];

					physicalglobals->K_row[index][tetrahedra[e].iIndex[j]] += (tmp);
					if (j > i){
						index = tetrahedra[e].iIndex[j];
						physicalglobals->K_row[index][tetrahedra[e].iIndex[i]] += (glm::transpose(tmp));
					}
				}
			}
			int idx = tetrahedra[e].iIndex[i];
			physicalglobals->F0[idx] -= Re * vTempForce;
		}
	}
}
//Clear the stiffness matrix
void clearStiffnessMatrix()
{
	for(size_t k=0; k<physicalglobals->total_points;k++){
		physicalglobals->F0[k].x = 0.0f;
		physicalglobals->F0[k].y = 0.0f;
		physicalglobals->F0[k].z = 0.0f;

		for(matrix_iterator Kij = physicalglobals->K_row[k].begin(); Kij != physicalglobals->K_row[k].end(); ++Kij){
			Kij->second = glm::mat3(0);
		}
	}
}
//Implement Gram_Schmidt orthogonalization
glm::mat3 Gram_Schmidt(glm::mat3 G)
{
	glm::vec3 row0(G[0][0], G[0][1], G[0][2]);
	glm::vec3 row1(G[1][0], G[1][1], G[1][2]);
	glm::vec3 row2(G[2][0], G[2][1], G[2][2]);

	float L0 = glm::length(row0);
	if(L0)
		row0 /= L0;
	row1 -= row0 * glm::dot(row0, row1);

	float L1 = glm::length(row1);
	if(L1)
		row1 /= L1;

	row2 = glm::cross( row0, row1);

	return glm::mat3(row0,
					 row1,
					 row2);
}
//Compute the orientation used for warping
void updateOrientation()
{
	for( size_t i=0; i < physicalglobals->iNumofTetrohedron; i++){
		float div6V = 1.0f / tetrahedra[i].fVolume * 6.0f;

		int i0 = tetrahedra[i].iIndex[0];
		int i1 = tetrahedra[i].iIndex[1];
		int i2 = tetrahedra[i].iIndex[2];
		int i3 = tetrahedra[i].iIndex[3];

		glm::vec3 p0 = physicalglobals->P[i0];
		glm::vec3 p1 = physicalglobals->P[i1];
		glm::vec3 p2 = physicalglobals->P[i2];
		glm::vec3 p3 = physicalglobals->P[i3];

		glm::vec3 e1 = p1 - p0;
		glm::vec3 e2 = p2 - p0;
		glm::vec3 e3 = p3 - p0;

		//G = [e1' e2' e3'][n1T n2T n3T]
		glm::vec3 n1 = glm::cross(e2,e3) * div6V;
		glm::vec3 n2 = glm::cross(e3,e1) * div6V;
		glm::vec3 n3 = glm::cross(e1,e2) * div6V;

		e1 = tetrahedra[i].e1;
		e2 = tetrahedra[i].e2;
		e3 = tetrahedra[i].e3;

		tetrahedra[i].Re[0][0] = e1.x * n1.x + e2.x * n2.x + e3.x * n3.x;
		tetrahedra[i].Re[0][1] = e1.x * n1.y + e2.x * n2.y + e3.x * n3.y;
		tetrahedra[i].Re[0][2] = e1.x * n1.z + e2.x * n2.z + e3.x * n3.z;

		tetrahedra[i].Re[1][0] = e1.y * n1.x + e2.y * n2.x + e3.y * n3.x;
		tetrahedra[i].Re[1][1] = e1.y * n1.y + e2.y * n2.y + e3.y * n3.y;
		tetrahedra[i].Re[1][2] = e1.y * n1.z + e2.y * n2.z + e3.y * n3.z;

		tetrahedra[i].Re[2][0] = e1.z * n1.x + e2.z * n2.x + e3.z * n3.x;
		tetrahedra[i].Re[2][1] = e1.z * n1.y + e2.z * n2.y + e3.z * n3.y;
		tetrahedra[i].Re[2][2] = e1.z * n1.z + e2.z * n2.z + e3.z * n3.z;

		tetrahedra[i].Re = Gram_Schmidt(tetrahedra[i].Re);
	}
}
//Reset the orientation to I
void resetOrientation()
{
	for(size_t i=0; i<physicalglobals->iNumofTetrohedron; i++){
		tetrahedra[i].Re = physicalglobals->eye;
	}
}
//Initialize the plasticity to 0
void initPlasticity()
{
	for(size_t i=0; i<tetrahedra.size(); i++){
		for(int j=0; j<6; j++){
			tetrahedra[i].fPlastic[j] = 0.0f;
		}
	}
}
//Compute Elasticity force
void forcePlasticity(float deltaT)
{
	for(size_t k=0; k<physicalglobals->iNumofTetrohedron;k++){
		float e_total[6];
		float e_elastic[6];
		for(unsigned int i=0; i<6; i++)
			e_elastic[i] = e_total[i] = 0.0f;

		//----------e_total = B((Re)^{-1}x-x0)
		for(unsigned int j=0; j<4; j++){
			glm::vec3 x_j = physicalglobals->P[tetrahedra[k].iIndex[j]];
			glm::vec3 x0_j = physicalglobals->Xi[tetrahedra[k].iIndex[j]];
			glm::mat3 ReT = glm::transpose(tetrahedra[k].Re);
			glm::vec3 prod = glm::vec3(ReT[0][0]*x_j.x + ReT[0][1]*x_j.y + ReT[0][2]*x_j.z,
									   ReT[1][0]*x_j.x + ReT[1][1]*x_j.y + ReT[1][2]*x_j.z,
									   ReT[2][0]*x_j.x + ReT[2][1]*x_j.y + ReT[2][2]*x_j.z);
			glm::vec3 tmp = prod - x0_j;

			float bj = tetrahedra[k].B[j].x;
			float cj = tetrahedra[k].B[j].y;
			float dj = tetrahedra[k].B[j].z;

			e_total[0] += bj*tmp.x;
			e_total[1] +=			 cj*tmp.y;
			e_total[2] +=					   dj*tmp.z;
			e_total[3] += cj*tmp.x + bj*tmp.y;
			e_total[4] += dj*tmp.x			  +bj*tmp.z;
			e_total[5] +=			 dj*tmp.y +cj*tmp.z;
		}

		//Elastic strain
		for(int i=0; i<6; i++)
			e_elastic[i] = e_total[i] - tetrahedra[k].fPlastic[i];

		//---------if elastic strain excedds cyield, then it is added to plastic strain by creep
		float fNormOfElastic = 0.0;
		for(int i=0; i<6; i++)
			fNormOfElastic += e_elastic[i]*e_elastic[i];
		fNormOfElastic = sqrt(fNormOfElastic);
		if(fNormOfElastic > yield){
			float amount = deltaT * MINIMUM(creep,1.0f / deltaT);
			for(int i=0; i<6; i++)
				tetrahedra[k].fPlastic[i] += amount*e_elastic[i];
		}

		//--------if plastic strain exceeds Cmax, it is clamped to maximum magnitude
		float fNormOfPlastic = 0.0;
		for(int i=0; i<6; i++)
			fNormOfPlastic += tetrahedra[k].fPlastic[i] * tetrahedra[k].fPlastic[i];
		fNormOfPlastic = sqrt(fNormOfPlastic);
		if(fNormOfPlastic > c_max){
			float fScale = c_max / fNormOfPlastic;
			for(int i=0; i<6; i++)
				tetrahedra[k].fPlastic[i] *= fScale;
		}
		for(size_t n=0; n<4; ++n){
			float* e_plastic = tetrahedra[k].fPlastic;

			float bn = tetrahedra[k].B[n].x;
			float cn = tetrahedra[k].B[n].y;
			float dn = tetrahedra[k].B[n].z;
			glm::vec3 f= glm::vec3(0);

			float bnD0 = bn*D0;
			float bnD1 = bn*D1;
			float bnD2 = bn*D2;
			float cnD0 = cn*D0;
			float cnD1 = cn*D1;
			float cnD2 = cn*D2;
			float dnD0 = dn*D0;
			float dnD1 = dn*D1;
			float dnD2 = dn*D2;

			f.x = bnD0*e_plastic[0] + bnD1*e_plastic[1] + bnD1*e_plastic[2] + cnD2*e_plastic[3] + dnD2*e_plastic[4];
			f.y = cnD1*e_plastic[0] + cnD0*e_plastic[1] + cnD1*e_plastic[2] + bnD2*e_plastic[3]					   + dnD2*e_plastic[5];
			f.z = dnD1*e_plastic[0] + dnD1*e_plastic[1] + dnD0*e_plastic[2] +					 bnD2*e_plastic[4] + cnD2*e_plastic[5];

			f *= tetrahedra[k].fVolume;
			int idx = tetrahedra[k].iIndex[n];
			physicalglobals->F[idx] += tetrahedra[k].Re*f;
		}
	}
}	
//Time integration
void dynamicsAssembly(float deltaT)
{
	float deltaT2 = deltaT * deltaT;
	for(size_t k=0; k<physicalglobals->total_points; k++){
		float m_i = physicalglobals->MASS[k];
		physicalglobals->b[k] = glm::vec3(0.0, 0.0, 0.0);

		matrix_map tmp = physicalglobals->K_row[k];
		matrix_iterator Kbegin = tmp.begin();
		matrix_iterator Kend = tmp.end();
		for(matrix_iterator K = Kbegin; K != Kend; K++){
			unsigned int j = K->first;
			glm::mat3 K_ij = K->second;
			glm::vec3 x_j = physicalglobals->P[j];
			glm::mat3& A_ij = physicalglobals->A_row[k][j];

			A_ij = K_ij * (deltaT2);
			glm::vec3 prod = glm::vec3(K_ij[0][0]*x_j.x + K_ij[0][1]*x_j.y + K_ij[0][2]*x_j.z,
									   K_ij[1][0]*x_j.x + K_ij[1][1]*x_j.y + K_ij[1][2]*x_j.z,
									   K_ij[2][0]*x_j.x + K_ij[2][1]*x_j.y + K_ij[2][2]*x_j.z);

			physicalglobals->b[k] -= prod;
			
			if(k == j){
				float c_i = fDamping*m_i;
				float tmp = m_i + deltaT*c_i;
				A_ij[0][0] += tmp;
				A_ij[1][1] += tmp;
				A_ij[2][2] += tmp;
			}
		}
		physicalglobals->b[k] -= physicalglobals->F0[k];
		physicalglobals->b[k] += physicalglobals->F[k];
		physicalglobals->b[k] *= deltaT;
		physicalglobals->b[k] += physicalglobals->V[k] * m_i;
	}
}
//Update positions
void updatePosition(float deltaT)
{
	for(size_t k=0; k<physicalglobals->total_points; k++){
		if(physicalglobals->IsFixed[k])
			continue;
		physicalglobals->P[k] += float(deltaT)*physicalglobals->V[k];
	}
}
//Collision detection
void collisionGround()
{
	for(size_t i =0; i<physicalglobals->total_points; i++){
		if(physicalglobals->P[i].y<0.0)
			physicalglobals->P[i].y=0.0;
	}
}
//Step Physics
void stepPhysics(float deltaT)
{
	computeforce();
	clearStiffnessMatrix();

	if(physicalglobals->bUsingStiffnessWarping)
		updateOrientation();
	else
		resetOrientation();

	stiffnessAssemble();
	forcePlasticity(deltaT);
	dynamicsAssembly(deltaT);
	conjugate_gradient_solver(deltaT);
	updatePosition(deltaT);
	//collisionGround();
}