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
	void halfStep(Sign sign);
	double integrate(Choice, unsigned int);

private:
	Representation q, p;
	double V(double);
	std::complex<double> *eV; // e^(-i/h V dt)
	std::complex<double> *eP; // e^(-i/2mh p^2 dt)
	double dt;
	const unsigned int size;
	const double invSize; // 1/N
	void setMomentum();
	fftw_plan pForward;
	fftw_plan pBackward;
};
