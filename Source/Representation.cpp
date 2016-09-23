#include "Representation.h"

Representation::Representation(unsigned int N)
{
	f = new std::complex<double>[N];
	x = new double[N];
}

Representation::~Representation()
{
	delete[] f;
	delete[] x;
}
