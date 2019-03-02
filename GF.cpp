#include "GF.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

CGF::CGF(void)
	: q(0)
	, p(0)
//	, PrimitivePoly(0)
//	, PrimitiveMatric(nullptr)
	, GFElement(nullptr)
	, TableAdd(nullptr)
	, TableMultiply(nullptr)
	, TableInverse(nullptr)
{
}


CGF::~CGF(void)
{
}


int CGF::GFAdd(int ele1, int ele2)
{
	return TableAdd[ele1][ele2];
}


int CGF::GFMultiply(int ele1, int ele2)
{
	return TableMultiply[ele1][ele2];
}


int CGF::GFInverse(int ele)
{
	if(ele == 0)
	{
		cerr << "Div 0 Error!" << endl;
		exit(-1);
	}
	return TableInverse[ele];
}


bool CGF::Initial(int GFq)
{
	// calculate order
	q = GFq;
	p = int(log(double(GFq)) / log(2.0));
	// allocate memory space
	GFElement = new CGFElement[q];
	for(int k = 0; k < q; k ++)
	{
		GFElement[k].ValueMatric = new int *[p]();
		for(int n = 0; n < p; n ++)
		{
			GFElement[k].ValueMatric[n] = new int[p]();
		}
	}
	TableAdd = new int *[q];
	for(int n = 0; n < q; n ++)
	{
		TableAdd[n] = new int [q];
	}
	TableMultiply = new int *[q];
	for(int n = 0; n < q; n ++)
	{
		TableMultiply[n] = new int [q];
	}
	TableInverse = new int [q];
	
	// read profile
	stringstream ss;
	ss << q << ".txt";
	//Arithmetic Table
	string ArithTableFileName = "./SRC/Arith.Table.GF.";
	ArithTableFileName += ss.str();
	ifstream ArithFin(ArithTableFileName);
	if (!ArithFin.is_open())
	{
		cerr << "Cannot open " << ArithTableFileName << endl ;
		exit(-1);
	}
	string rub;
	getline(ArithFin, rub);
//	cout << "Read Arithmetic Table File: " << rub << "..." << endl;
	ArithFin >> rub >> rub;
	for(int i = 0; i < q; i ++)
	{
		for(int j = 0; j < q; j ++)
		{
			ArithFin >> TableMultiply[i][j];
		}
	}
	ArithFin >> rub >> rub;
	for(int i = 0; i < q; i ++)
	{
		for(int j = 0; j < q; j ++)
		{
			ArithFin >> TableAdd[i][j];
		}
	}
	ArithFin >> rub >> rub;
	for(int i = 0; i < q; i ++)
	{
		ArithFin >> TableInverse[i];
	}
	ArithFin.close();
//	cout << "done." << endl;
	//Matric Representation
	string MatReprFileName = "./SRC/Mat.Repr.GF.";
	MatReprFileName += ss.str();
	ifstream MatReprFin(MatReprFileName);
	if (!MatReprFin.is_open())
	{
		cerr << "Cannot open " << MatReprFileName << endl;
		exit(-1);
	}
	getline(MatReprFin, rub);
//	cout << "Read Matric Representation File: " << rub << "..." << endl;
	// element zero is special issued
	GFElement[0].Order = q - 1;
	GFElement[0].ValuePoly = 0;
	for(int i = 0; i < p; i ++)
	{
		for(int j = 0; j < p; j ++)
		{
			GFElement[0].ValueMatric[i][j] = 0;
		}
	}
	// other non zero element
	for(int k = 0; k < q - 2; k ++)
	{
		int order, non0elepoly;
		MatReprFin >> rub >> rub >> rub >> order >> rub >> non0elepoly;
		GFElement[non0elepoly].Order = order;
		GFElement[non0elepoly].ValuePoly = non0elepoly;
		int tempInt2Char;
		for(int i = 0; i < p; i ++)
		{
			for(int j = 0; j < p; j++)
			{
				MatReprFin >> tempInt2Char;
				GFElement[non0elepoly].ValueMatric[i][j] = tempInt2Char;
			}
		}
	}
//	cout << "done." << endl;
	return true;
}
