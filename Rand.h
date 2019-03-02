#pragma once
#include <cmath>

class CRand
{
public:
	CRand(void);
	~CRand(void);
	unsigned long IX;
	unsigned long IY;
	unsigned long IZ;
	double Rand_Uniform(void);
	double Rand_Norm(double mu, double sigma);
};

