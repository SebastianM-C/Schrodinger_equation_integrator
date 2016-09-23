#include "Ket.h"
//#include "advisor-annotate.h"
#include <iostream>

Ket::Ket(double a, double b, double simTime, unsigned int N, unsigned int M)
		: q(N), p(N), size(N), invSize(1. / N)
{
	double dx = (b - a) / N;
	//	double dy = (b - a) / N;
	//	double dz = (b - a) / N;
	dt = simTime / M;

	pForward = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex *>(q.f), reinterpret_cast<fftw_complex *>(p.f), FFTW_FORWARD, FFTW_ESTIMATE);
	pBackward = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex *>(p.f), reinterpret_cast<fftw_complex *>(q.f), FFTW_BACKWARD, FFTW_ESTIMATE);

	//#pragma omp simd
	for (unsigned int i = 0; i < size; i++)
	{
		double x = a + i * dx;
		// double y = a + i * dx;
		// double z = a + i * dx;
		q.x[i] = x;

		if (i < N / 2)
		{
			p.x[i] = 2 * M_PI * i / (b - a);
		}
		else
		{
			p.x[i] = -2 * M_PI * (N - i) / (b - a);
		}

		// create initial probability distribution
		q.f[i] = 1 / (d * pow(M_PI, (double)0.25)) * exp(I * k0 * x - x * x / (2 * d * d));

		p.f[i] = 0;
	}
}

Ket::~Ket()
{
	fftw_destroy_plan(pForward);
	fftw_destroy_plan(pBackward);
}

// DEBUG
void Ket::print(std::ostream &out, Choice mode) const
{
	switch (mode)
	{
	case Q:
		for (unsigned int i = 0; i < size; i++)
		{
			out << q.x[i] << ' ' << std::norm(q.f[i]) << '\n';
		}
		break;
	case P:
		for (unsigned int i = 0; i < size; i++)
		{
			out << p.x[i] << ' ' << std::norm(p.f[i]) << '\n';
		}
		break;
	case Q_real_imag:
		for (unsigned int i = 0; i < size; i++)
		{
			out << q.x[i] << ' ' << std::real(q.f[i]) << ' ' << std::imag(q.f[i]) << '\n';
		}
		break;
	case P_real_imag:
		for (unsigned int i = 0; i < size; i++)
		{
			out << p.x[i] << ' ' << std::real(p.f[i]) << ' ' << std::imag(p.f[i]) << '\n';
		}
		break;
	}
}

void Ket::timeEvolution()
{
	//TODO: add potential energy operator

	fftw_execute(pForward);

	for (unsigned int i = 0; i < size; i++)
	{
		p.f[i] *= exp(-I * p.x[i] * p.x[i] * dt / (2. * m * h)); // apply time evolution operator
	}

	fftw_execute(pBackward);

	for (unsigned int i = 0; i < size; i++)
	{
		q.f[i] /= size; // normalize
	}
}

std::ostream &operator<<(std::ostream &out, const Ket &psi)
{
	unsigned int N = psi.size;
	for (unsigned int i = 0; i < N; i++)
	{
		out << psi.q.x[i] << ' ' << std::norm(psi.q.f[i]) << '\n';
	}
	return out;
}

void Ket::setMomentum()
{
	double delta = (q.x[size - 1] - q.x[0]) / size;
	for (unsigned int i = 0; i < size; i++)
	{
		p.f[i] /= exp(-I * p.x[i] * p.x[i] * dt / (2. * m * h)); // undo time evolution
		p.f[i] *= delta / sqrt(2. * M_PI * h);									 // normalize
	}
}

double Ket::norm(Choice choice)
{
	Representation *t;
	int first, last;
	if (choice == Q)
	{
		t = &q;
		first = 0;
		last = size - 1;
	}
	else
	{
		t = &p;
		first = size / 2;
		last = size / 2 - 1;
	}

	auto n = [&](int i) {
		return std::norm(t->f[i]);
	};
	double delta = (t->x[last] - t->x[first]) / size;
	double integral = delta / 2 * (n(first) + n(last));
	for (unsigned int i = 1; i < size - 1; i++)
	{
		integral += n(i) * delta;
	}

	return integral;
}

double Ket::mean(Choice choice)
{
	Representation *t;
	int first, last;
	if (choice == Q)
	{
		t = &q;
		first = 0;
		last = size - 1;
	}
	else
	{
		t = &p;
		first = size / 2;
		last = size / 2 - 1;
	}

	auto m = [&](int i) {
		return t->x[i] * std::norm(t->f[i]);
	};
	double delta = (t->x[last] - t->x[first]) / size;
	double integral = delta / 2 * (m(first) + m(last));
	for (unsigned int i = 1; i < size - 1; i++)
	{
		integral += m(i) * delta;
	}
	return integral;
}

double Ket::sqMean(Choice choice)
{
	Representation *t;
	int first, last;
	if (choice == Q)
	{
		t = &q;
		first = 0;
		last = size - 1;
	}
	else
	{
		t = &p;
		first = size / 2;
		last = size / 2 - 1;
	}

	auto sqm = [&](int i) {
		return t->x[i] * t->x[i] * std::norm(t->f[i]);
	};
	double delta = (t->x[last] - t->x[first]) / size;
	double integral = delta / 2 * (sqm(first) + sqm(last));
	for (unsigned int i = 1; i < size - 1; i++)
	{
		integral += sqm(i) * delta;
	}
	return integral;
}
