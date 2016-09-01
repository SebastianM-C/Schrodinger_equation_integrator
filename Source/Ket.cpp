#include "Ket.h"
//#include "advisor-annotate.h"

Ket::Ket(double a, double b, double simTime, unsigned int N, unsigned int M)
		: size(N)
{
	Representation *ri;
	Representation *pi;
	double dx = (b - a) / N;
	//	double dy = (b - a) / N;
	//	double dz = (b - a) / N;
	dt = simTime / M;

	r = new Representation *[N];
	p = new Representation *[N];

	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * size);
	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * size);

	pForward = fftw_plan_dft_1d(size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	pBackward = fftw_plan_dft_1d(size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

#pragma omp simd
	for (unsigned int i = 0; i < size; i++)
	{
		ri = new (std::nothrow) Representation;
		pi = new (std::nothrow) Representation;
		double x = a + i * dx;
		//double y = a + i*dx;
		//double z = a + i*dx;
		ri->x = x;
		// ri->y = y;
		// ri->z = z;

		if (i < N / 2)
		{
			pi->x = 2 * M_PI * i / (b - a);
		}
		else
		{
			pi->x = -2 * M_PI * (N - i) / (b - a);
		}
		// create initial probability distribution
		ri->fx = 1 / (d * pow(M_PI, (double)0.25)) * exp(I * k0 * x - x * x / (2 * d * d));

		r[i] = ri;
		p[i] = pi;
	}
}

Ket::~Ket()
{
	for (unsigned int i = 0; i < size; i++)
	{
		delete r[i];
	}
	for (unsigned int i = 0; i < size; i++)
	{
		delete p[i];
	}
	delete[] r;
	delete[] p;

	fftw_destroy_plan(pForward);
	fftw_destroy_plan(pBackward);
	fftw_free(in);
	fftw_free(out);
}

// DEBUG
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
	if (mode == 1)
	{
		for (unsigned int i = 0; i < size; i++)
		{
			out << r[i]->x << ' ' << r[i]->fx.real() << ' ' << r[i]->fx.imag() << '\n';
		}
	}
	if (mode == 2)
	{
		for (unsigned int i = 0; i < size; i++)
		{
			out << p[i]->x << ' ' << std::norm(p[i]->fx) << '\n';
		}
	}
}

void Ket::timeEvolution()
{
	for (unsigned int i = 0; i < size; i++)
	{
		memcpy(&in[i], &(r[i]->fx), sizeof(fftw_complex));
	}

	fftw_execute(pForward);

	for (unsigned int i = 0; i < size; i++)
	{
		memcpy(&(p[i]->fx), &out[i], sizeof(fftw_complex));
		p[i]->fx *= exp(-I * p[i]->x * p[i]->x * dt / (2. * m * h));
		memcpy(&in[i], &(p[i]->fx), sizeof(fftw_complex));
	}

	fftw_execute(pBackward);

	for (unsigned int i = 0; i < size; i++)
	{
		memcpy(&(r[i]->fx), &out[i], sizeof(fftw_complex));
		r[i]->fx /= size;
	}
}

std::ostream &operator<<(std::ostream &out, const Ket &psi)
{
	unsigned int N = psi.size;
	for (unsigned int i = 0; i < N; i++)
	{
		out << psi.r[i]->x << ' ' << std::norm(psi.r[i]->fx) << '\n';
	}
	return out;
}

double Ket::mean()
{
	auto m = [&](int i) {
		return r[i]->x * std::norm(r[i]->fx);
		//return p[i]->x * std::norm(p[i]->fx);
	};
	double delta = (r[size - 1]->x - r[0]->x) / size;
	//double delta = (p[size - 1]->x - p[0]->x) / size;
	double integral = delta / 2 * (m(0) + m(size - 1));
	for (unsigned int i = 1; i < size - 1; i++)
	{
		integral += m(i) * delta;
	}
	return integral;
}

double Ket::sqMean()
{
	auto sqm = [&](int i) {
		return r[i]->x * r[i]->x * std::norm(r[i]->fx);
		//return p[i]->x * p[i]->x * std::norm(p[i]->fx);
	};
	double delta = (r[size - 1]->x - r[0]->x) / size;
	//double delta = (p[size - 1]->x - p[0]->x) / size;
	double integral = delta / 2 * (sqm(0) + sqm(size - 1));
	for (unsigned int i = 1; i < size - 1; i++)
	{
		integral += sqm(i) * delta;
	}
	return integral;
}
