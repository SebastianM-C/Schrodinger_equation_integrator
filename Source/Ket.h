#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <ostream>
#include "position.h"
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

private:
	Position **r;
	double dt;
	const unsigned int size;
};
