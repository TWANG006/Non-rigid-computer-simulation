#ifndef SPRING_H
#define SPRING_H

class Spring
{
public:
	int iP1,iP2;
	float fRestlength;
	float fKs,fKd;		//Spring stiffness and damping coefficients.
	int iType;
};

#endif // !SPRING_H

