#include "Ket.h"
//#include "advisor-annotate.h"
#include <iostream>

Ket::Ket(double a, double b, double simTime, unsigned int N, unsigned int M)
		: size(N)
{
	Representation *qi;
	Representation *pi;
	double dx = (b - a) / N;
	//	double dy = (b - a) / N;
	//	double dz = (b - a) / N;
	dt = simTime / M;

	q = new Representation *[N];
	p = new Representation *[N];

	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * size);
	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * size);

	pForward = fftw_plan_dft_1d(size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	pBackward = fftw_plan_dft_1d(size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

	//#pragma omp simd
	for (unsigned int i = 0; i < size; i++)
	{
		qi = new (std::nothrow) Representation;
		pi = new (std::nothrow) Representation;
		double x = a + i * dx;
		//double y = a + i*dx;
		//double z = a + i*dx;
		qi->x = x;
		// qi->y = y;
		// qi->z = z;

		if (i < N / 2)
		{
			pi->x = 2 * M_PI * i / (b - a);
		}
		else
		{
			pi->x = -2 * M_PI * (N - i) / (b - a);
		}
		// create initial probability distribution
		qi->fx = 1 / (d * pow(M_PI, (double)0.25)) * exp(I * k0 * x - x * x / (2 * d * d));

		q[i] = qi;
		p[i] = pi;
	}
}

Ket::~Ket()
{
	for (unsigned int i = 0; i < size; i++)
	{
		delete q[i];
	}
	for (unsigned int i = 0; i < size; i++)
	{
		delete p[i];
	}
	delete[] q;
	delete[] p;

	fftw_destroy_plan(pForward);
	fftw_destroy_plan(pBackward);
	fftw_free(in);
	fftw_free(out);
}

// DEBUG
void Ket::print(std::ostream &out, Choice mode) const
{
	switch (mode)
	{
	case Q:
		for (unsigned int i = 0; i < size; i++)
		{
			out << q[i]->x << ' ' << std::norm(q[i]->fx) << '\n';
		}
		break;
	case P:
		for (unsigned int i = 0; i < size; i++)
		{
			out << p[i]->x << ' ' << std::norm(p[i]->fx) << '\n';
		}
		break;
	case Q_real_imag:
		for (unsigned int i = 0; i < size; i++)
		{
			out << q[i]->x << ' ' << std::real(q[i]->fx) << ' ' << std::imag(q[i]->fx) << '\n';
		}
		break;
	case P_real_imag:
		for (unsigned int i = 0; i < size; i++)
		{
			out << p[i]->x << ' ' << std::real(p[i]->fx) << ' ' << std::imag(p[i]->fx) << '\n';
		}
		break;
	}
}

void Ket::timeEvolution()
{
	for (unsigned int i = 0; i < size; i++)
	{
		memcpy(&in[i], &(q[i]->fx), sizeof(fftw_complex));
	}

	fftw_execute(pForward);

	for (unsigned int i = 0; i < size; i++)
	{
		memcpy(&(p[i]->fx), &out[i], sizeof(fftw_complex));
		p[i]->fx *= exp(-I * p[i]->x * p[i]->x * dt / (2. * m * h)); // apply time evolution operator
		memcpy(&in[i], &(p[i]->fx), sizeof(fftw_complex));
	}

	fftw_execute(pBackward);

	for (unsigned int i = 0; i < size; i++)
	{
		memcpy(&(q[i]->fx), &out[i], sizeof(fftw_complex));
		q[i]->fx /= size; // normalize
	}
}

std::ostream &operator<<(std::ostream &out, const Ket &psi)
{
	unsigned int N = psi.size;
	for (unsigned int i = 0; i < N; i++)
	{
		out << psi.q[i]->x << ' ' << std::norm(psi.q[i]->fx) << '\n';
	}
	return out;
}

void Ket::setMomentum()
{
	double delta = (q[size - 1]->x - q[0]->x) / size;
	//std::cout<<delta<<'\n';
	for (unsigned int i = 0; i < size; i++)
	{
		p[i]->fx /= exp(-I * p[i]->x * p[i]->x * dt / (2. * m * h)); // undo time evolution
		p[i]->fx *= delta / sqrt(2. * M_PI * h);										 // normalize
	}
}

double Ket::norm(Choice choice)
{
	Representation **t;
	int first, last;
	if (choice == Q)
	{
		t = q;
		first = 0;
		last = size - 1;
	}
	else
	{
		t = p;
		first = size / 2;
		last = size / 2 - 1;
	}

	auto n = [&](int i)
	{
		return std::norm(t[i]->fx);
	};
	double delta = (t[last]->x - t[first]->x) / size;
	double integral = delta / 2 * (n(first) + n(last));
	for (unsigned int i = 1; i < size - 1; i++)
	{
		integral += n(i) * delta;
	}

	return integral;
}

double Ket::mean(Choice choice)
{
	Representation **t;
	int first, last;
	if (choice == Q)
	{
		t = q;
		first = 0;
		last = size - 1;
	}
	else
	{
		t = p;
		first = size / 2;
		last = size / 2 - 1;
	}

	auto m = [&](int i)
	{
		return t[i]->x * std::norm(t[i]->fx);
	};
	double delta = (t[last]->x - t[first]->x) / size;
	double integral = delta / 2 * (m(first) + m(last));
	for (unsigned int i = 1; i < size - 1; i++)
	{
		integral += m(i) * delta;
	}
	return integral;
}

double Ket::sqMean(Choice choice)
{
	Representation **t;
	int first, last;
	if (choice == Q)
	{
		t = q;
		first = 0;
		last = size - 1;
	}
	else
	{
		t = p;
		first = size / 2;
		last = size / 2 - 1;
	}

	auto sqm = [&](int i)
	{
		return t[i]->x * t[i]->x * std::norm(t[i]->fx);
	};
	double delta = (t[last]->x - t[first]->x) / size;
	double integral = delta / 2 * (sqm(first) + sqm(last));
	for (unsigned int i = 1; i < size - 1; i++)
	{
		integral += sqm(i) * delta;
	}
	return integral;
}
