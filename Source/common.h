#pragma once

#include <complex>

// constants
extern const double h;
extern const double d;
extern const double k0;
extern const double m;
extern const double inv2mh;
extern const std::complex<double> I;

enum Choice
{
	Q,
	P,
	Q_real_imag,
	P_real_imag
};
