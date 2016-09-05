/**
 * Schrodinger equation integrator
 *
 */

#include <fstream>
#include <iostream>

#include "common.h"
#include "ket.h"

// constants
const double h = 1;
const double d = 1;
const double k0 = 1;
const double m = 1;
const std::complex<double> I(0., 1.);

int main(int argc, char *argv[])
{
	std::ifstream in("input.dat");
	std::ofstream out("output.dat");
	double a, b, simTime, t = 0;
	unsigned int N, M;
	in >> a >> b >> simTime >> N >> M;
	Ket Psi(a, b, simTime, N, M);

	auto variance = [&Psi](bool flag)
	{
		return Psi.sqMean(flag) - (Psi.mean(flag) * Psi.mean(flag));
	};
	auto uncertainty = [&Psi, &variance]()
	{
		return variance(1) * variance(0);
	};
	//Psi.print(out, 0);
	while (t < simTime)
	{
		Psi.timeEvolution();
		t += (simTime / M);
		Psi.setMomentum();
		out << t << ' ' << variance(0) << '\n';
		//out << t << ' ' << uncertainty() << '\n';
		//Psi.print(out, 1);
		//out << "\n\n";
	}

	return 0;
}
