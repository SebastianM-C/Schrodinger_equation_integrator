#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <ostream>
#include <cstring>
#include "Representation.h"
#include "common.h"
#include <fftw3.h>

class Ket
{
public:
	Ket(double, double, double, unsigned int, unsigned int);
	~Ket();
	void print(std::ostream &, int = 0) const;
	friend std::ostream &operator<<(std::ostream &, const Ket &);
	void timeEvolution();
	void setMomentum();
	double mean(bool);
	double sqMean(bool);

private:
	Representation **r, **p;
	double dt;
	const unsigned int size;
	fftw_plan pForward;
	fftw_plan pBackward;
	fftw_complex *in, *out;
};
