#include "Rand.h"


CRand::CRand(void)
	: IX(0)
	, IY(0)
	, IZ(0)
{
}


CRand::~CRand(void)
{
}


double CRand::Rand_Uniform(void)
{
	double temp = 0.0;
	IX = (IX * 249) % 61967;
	IY = (IY * 251) % 63443;
	IZ = (IZ * 252) % 63599;
	temp = (((double)IX) / ((double)61967)) + (((double)IY) / ((double)63443))
		+ (((double)IZ) / ((double)63599));
	temp -= (int)temp;
	//double temp = double(rand() % 32767) / 32767.0;
	return temp;
}


double CRand::Rand_Norm(double mu, double sigma)
{
	double u1, u2, u;
	u1 = Rand_Uniform();
	u2 = Rand_Uniform();
	u = mu + sigma* cos(2 * acos(-1.0) * u2) * sqrt(-2.0 * log(1.0 - u1));
	return u;
}
