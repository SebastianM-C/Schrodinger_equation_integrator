#pragma once

#include <complex>

struct Representation
{
	Representation(unsigned int);
	~Representation();
	std::complex<double> *f; // function values
	double *x;							 // coordinates
};
