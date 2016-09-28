#include "Ket.h"
#include <iostream>

Ket::Ket(double a, double b, double simTime, unsigned int N, unsigned int M)
		: q(N), p(N), size(N), invSize(1. / N)
{
	double dx = (b - a) / N;
	//	double dy = (b - a) / N;
	//	double dz = (b - a) / N;
	dt = simTime / M;

	eV = new std::complex<double>[N]; // e^(-i/h V dt)
	eP = new std::complex<double>[N]; // e^(-i/2mh p^2 dt)

	pForward = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex *>(q.f),
															reinterpret_cast<fftw_complex *>(p.f), FFTW_FORWARD, FFTW_ESTIMATE);
	pBackward = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex *>(p.f),
															 reinterpret_cast<fftw_complex *>(q.f), FFTW_BACKWARD, FFTW_ESTIMATE);

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

		// potential energy operator values
		eV[i] = exp(-I / h * V(x) * dt);
		// kinetic energy operator values
		eP[i] = exp(-I * p.x[i] * p.x[i] * dt * inv2mh);
	}
}

Ket::~Ket()
{
	delete[] eV;
	delete[] eP;
	fftw_destroy_plan(pForward);
	fftw_destroy_plan(pBackward);
}

// potential energy
double Ket::V(double x)
{
	double x0 = 2, V0 = 1.5;
	if (x < x0)
		return 0;
	return V0;
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

// Integrate the Schrodinger equation using the split operator method
void Ket::timeEvolution()
{
	fftw_execute(pForward); // Fourier transform into momentum space
#pragma omp simd
	for (unsigned int i = 0; i < size; i++)
	{
		p.f[i] *= eP[i]; // apply kinetic energy operator
	}

	fftw_execute(pBackward); // Fourier transform back into position space
#pragma omp simd
	for (unsigned int i = 0; i < size; i++)
	{
		q.f[i] *= eV[i] * invSize; // apply potential energy operator and normalize
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

void Ket::halfStep(Sign sign)
{ // Forward = 0, Backward = 2
	//#pragma omp simd
	for (unsigned int i = 0; i < size; i++)
	{ // advance half step for -1 or go back for 1
		std::complex<double> ev = exp(((double)sign - 1) * I / h * V(q.x[i]) * 0.5 * dt);
		q.f[i] *= ev;
	}
}

void Ket::setMomentum()
{
	double delta = (q.x[size - 1] - q.x[0]) / size;
#pragma omp simd
	for (unsigned int i = 0; i < size; i++)
	{
		p.f[i] *= exp(I * p.x[i] * p.x[i] * dt * inv2mh); // undo time evolution
		p.f[i] *= delta * invsqrt2pih;										// normalize
	}
}

// compute integral x^power * choice(x) dx
double Ket::integrate(Choice choice, unsigned int power)
{
	Representation *t;
	int first, last;
	halfStep(Backward);
	if (choice == Q)
	{
		t = &q;
		first = 0;
		last = size - 1;
	}
	else
	{
		setMomentum();
		t = &p;
		first = size / 2;
		last = size / 2 - 1;
	}

	auto f = [&](int i)
	{
		double x = 1;
		for (unsigned int j = 0; j < power; j++)
		{
			x *= t->x[i];
		}
		return x * std::norm(t->f[i]);
	};
	double delta = (t->x[last] - t->x[first]) / size;
	double integral = delta / 2 * (f(first) + f(last));
#pragma omp simd
	for (unsigned int i = 1; i < size - 1; i++)
	{
		integral += f(i) * delta;
	}

	halfStep(Forward);
	return integral;
}
