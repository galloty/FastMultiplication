/*
Copyright 2024, Yves Gallot

FastMultiplication is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <complex>
#include <iostream>

typedef std::complex<double> Complex;

// P = P^2 mod R, where R = x^m - r
static void square_fast(Complex * const P, const size_t m, const Complex & r)
{
	if (m == 1) { P[0] *= P[0]; return; } 	// P is a scalar

	// R = R0 * R1, where Ri = x^{m/2} - ri
	const Complex r0 = std::sqrt(r);	// r1 = -r0;

	// Pi = P mod Ri: if P = a[0] + a[1].x + ... then Pi mod (x^{m/2} - ri) = (a[0] + ri * a[m/2]) + ...
	Complex * const P0 = &P[0 * m / 2]; Complex * const P1 = &P[1 * m / 2];

	for (size_t i = 0; i < m / 2; ++i)
	{
		const Complex u0 = P0[i], u1 = r0 * P1[i];
		P0[i] = u0 + u1;
		P1[i] = u0 - u1;
	}

	// P1 = P^2 mod R1, P2 = P^2 mod R2 => divide and conquer
	square_fast(P0, m / 2, r0); square_fast(P1, m / 2, -r0);

	// Chinese remainder theorem or simply the inverse of direct butterfly
	for (size_t i = 0; i < m / 2; ++i)
	{
		const Complex u0 = P0[i], u1 = P1[i];
		P0[i] = (u1 + u0) /  2.0;
		P1[i] = (u0 - u1) / (2.0 * r0);
	}
}

int main()
{
	const size_t n = 8;
	Complex P[n]; for (size_t i = 0; i < n; ++i) P[i] = (i + 1) * (i + 3);
	Complex r = 7;	// modulo x^n - r

	std::cout << "Mod(80*x^7 + 63*x^6 + 48*x^5 + 35*x^4 + 24*x^3 + 15*x^2 + 8*x + 3, x^8 - 7)^2" << std::endl;
	std::cout << " = 4608*x^7 + 47572*x^6 + 72128*x^5 + 82362*x^4 + 81920*x^3 + 74032*x^2 + 61536*x + 46902" << std::endl;

	square_fast(P, n, r);

	std::cout << " = ";
	for (size_t i = n - 1, j = 0; j < n; --i, ++j)
	{
		std::cout << P[i].real();
		if (i > 0)
		{
			std::cout << "*x";
			if (i > 1) std::cout << "^" << i;
			std::cout << " + ";
		}
	}
	std::cout << std::endl;

	return EXIT_SUCCESS;
}
