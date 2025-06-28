/*
Copyright 2024, Yves Gallot

FastMultiplication is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <iostream>

typedef std::complex<double> Complex;

// P = P^2 mod R, where R = x^m - r, r = e^{2*i*Pi * j * m / n}
static void square(Complex * const P, const size_t j, const size_t m, const size_t n)
{
	if (m == 1) { P[0] *= P[0]; return; }	// P is a scalar

	// R = R1 * R2, where R1 = x^{m/2} - sqrt_r, R2 = x^{m/2} + sqrt_r
	const Complex sqrt_r = std::polar(1.0, M_PI * j * m / n);

	// P1 = P mod R1, P2 = P mod R2
	// If P = a[0] + a[1].x + ... then P mod (x^{m/2} - sr) = (a[0] + a[m/2] * sr) + ...
	Complex * const P1 = &P[0 * m / 2]; Complex * const P2 = &P[1 * m / 2];
	for (size_t i = 0; i < m / 2; ++i)
	{
		const Complex u0 = P1[i], u1 = P2[i];
		P1[i] = u0 + sqrt_r * u1;
		P2[i] = u0 - sqrt_r * u1;
	}
	
	// P1 = P^2 mod R1, P2 = P^2 mod R2 => divide and conquer
	square(P1, j, m / 2, n); square(P2, j + n / m, m / 2, n);	// If j = n / m then sqrt_r = -1

	// Chinese remainder theorem or simply the inverse of direct butterflies
	for (size_t i = 0; i < m / 2; ++i)
	{
		const Complex u0 = P1[i], u1 = P2[i];
		P1[i] = 0.5 * (u0 + u1);
		P2[i] = 0.5 * (u0 - u1) * std::conj(sqrt_r);	// The inverse of sqrt_r is its conjugate
	}
}

// P = P^2 mod R, where R = x^m - 1
static void square_twisted(Complex * const P, const size_t m)	// j = 0, n is not needed
{
	if (m == 1) { P[0] *= P[0]; return; }

	// The root is one
	Complex * const P1 = &P[0 * m / 2]; Complex * const P2 = &P[1 * m / 2];
	for (size_t i = 0; i < m / 2; ++i)
	{
		const Complex u0 = P1[i], u1 = P2[i];
		P1[i] = u0 + u1;
		P2[i] = u0 - u1;
	}

	// R2 = x^{m/2} + 1, P2 must be twisted such that R2' = x^{m/2} - 1
	for (size_t i = 0; i < m / 2; ++i) P2[i] *= std::polar(1.0, M_PI * i / (m / 2));

	square_twisted(P1, m / 2); square_twisted(P2, m / 2);

	for (size_t i = 0; i < m / 2; ++i)
	{
		const Complex u0 = P1[i], u1 = P2[i] * std::polar(1.0, -M_PI * i / (m / 2));	// untwist
		P1[i] = 0.5 * (u0 + u1);
		P2[i] = 0.5 * (u0 - u1);
	}
}

static void display(const size_t n, const Complex * const P)
{
	std::cout << "       ";
	for (size_t i = n - 1, j = 0; j < n; --i, ++j)
	{
		std::cout << P[i].real();
		if (i > 0)
		{
			std::cout << " x";
			if (i > 1) std::cout << "^" << i;
			std::cout << " + ";
		}
	} 
	std::cout << std::endl;
}

int main()
{
	const size_t n = 16;
	Complex P[n];

	std::cout << "Mod(1 + 2*x + 3*x^2 + 6*x^3 + 5*x^4 + 10*x^5 + 7*x^6 + 14*x^7 + 9*x^8 + 18*x^9 + 11*x^10 + 22*x^11 + 13*x^12 + 26*x^13 + 15*x^14 + 30*x^15, x^16 - 1)^2" << std::endl;
	std::cout << " = Mod(1376*x^15 + 2168*x^14 + 1824*x^13 + 2600*x^12 + 2144*x^11 + 2872*x^10 + 2336*x^9 + 2984*x^8 + 2400*x^7 + 2936*x^6 + 2336*x^5 + 2728*x^4 + 2144*x^3 + 2360*x^2 + 1824*x + 1832, x^16 - 1)" << std::endl;

	for (size_t i = 0; i < n; ++i) P[i] = double((i % 2 == 0) ? i + 1 : i * 2);
	// modulo x^n - 1
	square(P, 0, n, n);	// 1 = e^{2*i*Pi * 0 / n}
	display(n, P);

	for (size_t i = 0; i < n; ++i) P[i] = double((i % 2 == 0) ? i + 1 : i * 2);
	square_twisted(P, n);
	display(n, P);

	return EXIT_SUCCESS;
}
