#include "conjugate_gradient.h"
#include "physicalSimulation.h"

/*void initialize_CGsolver()
{
	cgglobals = new ConjugateGradient();

}*/

ConjugateGradient* cgglobals;

//Conjugate gradient solver
void conjugate_gradient_solver(float dt)
{
	int i=0;
	for(size_t k=0; k < physicalglobals->total_points; k++){
		if(physicalglobals->IsFixed[k])
			continue;
		cgglobals->residual[k] = physicalglobals->b[k];

		matrix_iterator Abegin = physicalglobals->A_row[k].begin();
		matrix_iterator Aend = physicalglobals->A_row[k].end();
		for(matrix_iterator it_A = Abegin; it_A != Aend; ++it_A){
			unsigned int j= it_A->first;
			glm::mat3& A_ij = it_A->second;
			float v_jx = physicalglobals->V[j].x;
			float v_jy = physicalglobals->V[j].y;
			float v_jz = physicalglobals->V[j].z;
			glm::vec3 prod = glm::vec3(A_ij[0][0]*v_jx + A_ij[0][1]*v_jy + A_ij[0][2]*v_jz,
									   A_ij[1][0]*v_jx + A_ij[1][1]*v_jy + A_ij[1][2]*v_jz,
									   A_ij[2][0]*v_jx + A_ij[2][1]*v_jy + A_ij[2][2]*v_jz);

			cgglobals->residual[k] -= prod;
		}
		cgglobals->d[k] = cgglobals->residual[k];
	}
	for(int i=0; i<i_max; i++){
		float d0 = 0;
		float dTq = 0;
		for(size_t k=0; k<physicalglobals->total_points;k++){
			if(physicalglobals->IsFixed[k])
				continue;
			cgglobals->q[k].x = 0.0; cgglobals->q[k].y = 0.0; cgglobals->q[k].z = 0.0;

			matrix_iterator Abegin = physicalglobals->A_row[k].begin();
			matrix_iterator Aend   = physicalglobals->A_row[k].end();
			for(matrix_iterator it_A = Abegin; it_A != Aend; ++it_A){
				unsigned int j = it_A->first;
				glm::mat3& A_ij = it_A->second;
				float dX = cgglobals->d[j].x;
				float dY = cgglobals->d[j].y;
				float dZ = cgglobals->d[j].z;
				glm::vec3 prod = glm::vec3(A_ij[0][0]*dX + A_ij[0][1]*dY + A_ij[0][2]*dZ,
										   A_ij[1][0]*dX + A_ij[1][1]*dY + A_ij[1][2]*dZ,
										   A_ij[2][0]*dX + A_ij[2][1]*dY + A_ij[2][2]*dZ);
				cgglobals->q[k] += prod;
			}
			d0 += glm::dot(cgglobals->residual[k],cgglobals->residual[k]);
			dTq  += glm::dot(cgglobals->d[k], cgglobals->q[k]);
		}
		if(fabs(dTq) < tiny)
			dTq = tiny;

		float alfa = d0 / dTq;
		float dnew = 0.0;

		for(size_t k=0; k<physicalglobals->total_points; k++){
			if(physicalglobals->IsFixed[k])
				continue;

			physicalglobals->V[k] += cgglobals->d[k]*alfa;
			cgglobals->residual[k] -= cgglobals->q[k]*alfa;
			dnew += glm::dot(cgglobals->residual[k], cgglobals->residual[k]);
		}
		if( i >= i_max && dnew <tolerence)
			break;
		if(fabs(d0) < tiny)
			d0 = tiny;

		float beta = dnew / d0;

		for(size_t k=0; k<physicalglobals->total_points; k++){
			if(physicalglobals->IsFixed[k])
				continue;
			cgglobals->d[k] = cgglobals->residual[k] + cgglobals->d[k]*beta;
		}
	}
}