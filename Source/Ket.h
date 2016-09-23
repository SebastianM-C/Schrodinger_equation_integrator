#pragma once

#define _USE_MATH_DEFINES
#include "Representation.h"
#include "common.h"
#include <cmath>
#include <cstring>
#include <fftw3.h>
#include <ostream>

class Ket
{
public:
	Ket(double, double, double, unsigned int, unsigned int);
	~Ket();
	void print(std::ostream &, Choice = Q) const;
	friend std::ostream &operator<<(std::ostream &, const Ket &);
	void timeEvolution();
	void setMomentum();
	double norm(Choice);
	double mean(Choice);
	double sqMean(Choice);

private:
	Representation q, p;
	double dt;
	const unsigned int size;
	const double invSize;			// 1/N
	//std::complex<double> *eV; // e^(V)
	fftw_plan pForward;
	fftw_plan pBackward;
};
