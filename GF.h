#pragma once
#include "GFElement.h"

class CGF
{
public:
	CGF(void);
	~CGF(void);
	int q;
	int p;
//	int PrimitivePoly;
//	char** PrimitiveMatric;
	CGFElement* GFElement;
	int** TableAdd;
	int** TableMultiply;
	int* TableInverse;
public:
	int GFAdd(int ele1, int ele2);
	int GFMultiply(int ele1, int ele2);
	int GFInverse(int ele);
	bool Initial(int GFq);
};

