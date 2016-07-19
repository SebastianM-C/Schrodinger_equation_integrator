#include "ket.h"
//#include "advisor-annotate.h"

Ket::Ket(double a, double b, double simTime, unsigned int N, unsigned int M)
:size(N)
{
	Position *p;
	double dx = (b - a) / N;
//	double dy = (b - a) / N;
//	double dz = (b - a) / N;
	dt = simTime / M;
	r = new Position*[N];

	#pragma omp simd
	for (unsigned int i = 0; i < size; i++)
	{
		p = new (std::nothrow) Position;
		double x = a + i*dx;
		//double y = a + i*dx;
		//double z = a + i*dx;
		p->x = x;
		//p->y = y;
		//p->z = z;

		// create initial probability distribution
		p->fx.real(1 / (d*pow(M_PI, (double)0.25)) * exp(-x*x / (2 * d*d)) * cos(k0*x));
		p->fx.imag(1 / (d*pow(M_PI, (double)0.25)) * exp(-x*x / (2 * d*d)) * sin(k0*x));
		
		r[i] = p;
	}
}

Ket::~Ket()
{
	for (unsigned int i = 0; i < size; i++)
	{
		delete r[i];
	}
	delete[] r;
}

void Ket::print(std::ostream &out, int mode) const
{
	if (mode == 0)
	{
		for (unsigned int i = 0; i < size; i++)
		{
			out << r[i]->x << ' ' << std::norm(r[i]->fx) << '\n';
		}
		return;
	}
	//if (mode == 1)
	{
		for (unsigned int i = 0; i < size; i++)
		{
			out << r[i]->x << ' ' << r[i]->fx.real() << ' ' << r[i]->fx.imag() << '\n';
		}
	}
}

void Ket::timeEvolution()
{
	fftw_plan p;
	fftw_complex *in, *out;
	int N = size;

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	fftw_execute(p);

	for (unsigned int i = 0; i < N; i++)
	{
		// TODO: something 
	}

	fftw_destroy_plan(p);
	fftw_free(in); fftw_free(out);
}



std::ostream& operator << (std::ostream& out, const Ket& psi)
{
	int N = psi.size;
	for (unsigned int i = 0; i < N; i++)
	{
		out << psi.r[i]->x << ' ' << std::norm(psi.r[i]->fx) << '\n';
	}
	return out;
}
