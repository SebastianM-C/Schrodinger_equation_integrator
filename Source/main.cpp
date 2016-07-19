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

int main(int argc, char *argv[])
{
	std::ifstream in("input.dat");
	std::ofstream out("output.dat");
	double a, b, simTime;
	unsigned int N, M;
	in >> a >> b >> simTime >> N >> M;
	Ket Psi(a, b, simTime, N, M);

	for (unsigned int i = 0; i < M; i++)
	{
		Psi.timeEvolution();
		//out << Psi;
		Psi.print(out);
		out << '\n';
		break;
	}
	
	return 0;
}
