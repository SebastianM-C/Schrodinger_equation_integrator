/**
 * Schrodinger equation integrator
 *
 */

#include "Ket.h"
#include "common.h"
#include <fstream>
#include <iostream>

// constants
const double h = 1;
const double d = 1;
const double k0 = 1;
const double m = 1;
const double inv2mh = 1. / (2. * m * h);
const double invsqrt2pih = 1. / sqrt(2. * M_PI * h);
const std::complex<double> I(0., 1.);

int main(int argc, char *argv[])
{
	std::ifstream in("input.dat");
	std::ofstream out("output.dat");
	double a, b, simTime, dt, t = 0;
	unsigned int N, M;
	in >> a >> b >> simTime >> N >> M;
	Ket Psi(a, b, simTime, N, M);
	dt = simTime / M;

	auto variance = [&Psi](Choice flag)
	{
		return Psi.sqMean(flag) - (Psi.mean(flag) * Psi.mean(flag));
	};

	Psi.setPotential();
	while (t < simTime)
	{
		Psi.timeEvolution();
		t += dt;
		out << t << ' ' << Psi.norm(Q) << '\n';
		//out << t << ' ' << variance(P) * variance(Q) << '\n';
		//Psi.print(out, Q);
		//out << "\n\n";
	}
	Psi.setPotential();

	return 0;
}
