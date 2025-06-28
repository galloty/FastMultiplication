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
static void square(Complex * const P, const Complex & r, const size_t m)
{
	if (m == 1) { P[0] *= P[0]; return; } 	// P is a scalar

	// R = R1 * R2, where R1 = x^{m/2} - sr, R2 = x^{m/2} + sr
	const Complex sr = std::sqrt(r);

	// P1 = P mod R1, P2 = P mod R2
	// If P = a[0] + a[1].x + ... then P mod (x^{m/2} - sr) = (a[0] + a[m/2] * sr) + ...
	Complex * const P1 = &P[0 * m/2]; Complex * const P2 = &P[1 * m/2];
	for (size_t i = 0; i < m/2; ++i)
	{
		const Complex u0 = P1[i], u1 = P2[i];
		P1[i] = u0 + sr * u1;
		P2[i] = u0 - sr * u1;
	}
	
	// P1 = P^2 mod R1, P2 = P^2 mod R2 => divide and conquer
	square(P1, sr, m/2); square(P2, -sr, m/2);

	// Chinese remainder theorem or simply the inverse of direct butterfly
	const Complex sr_inv = 1.0 / sr;
	for (size_t i = 0; i < m/2; ++i)
	{
		const Complex u0 = P1[i], u1 = P2[i];
		P1[i] = 0.5 * (u0 + u1);
		P2[i] = 0.5 * (u0 - u1) * sr_inv;
	}
}

int main()
{
	const size_t n = 8;
	// P = 1 + 2*x + 3*x^2 + 6*x^3 + 5*x^4 + 10*x^5 + 7*x^6 + 14*x^7;
	Complex P[n]; for (size_t i = 0; i < n; ++i) P[i] = double((i % 2 == 0) ? i + 1 : i * 2);
	// modulo x^n - 1
	Complex r = 2;
	square(P, r, n);

	std::cout << "Mod(1 + 2*x + 3*x^2 + 6*x^3 + 5*x^4 + 10*x^5 + 7*x^6 + 14*x^7, x^8 - 2)^2" << std::endl;
	std::cout << " = Mod(176*x^7 + 512*x^6 + 468*x^5 + 701*x^4 + 584*x^3 + 686*x^2 + 540*x + 487, x^8 - 2)" << std::endl;
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

	return EXIT_SUCCESS;
}
