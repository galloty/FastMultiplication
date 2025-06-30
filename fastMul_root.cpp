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
#include <cstdlib>
#include <ctime>

typedef std::complex<double> Complex;

// P = P^2 mod R, where R = x^m - r
static void square_fast(Complex * const P, const size_t m, const Complex & r)
{
	if (m == 1) { P[0] *= P[0]; return; } 	// P is a scalar

	if (m % 3 == 0)
	{
		// R = R1 * R2 * r3, where Ri = x^{m/3} - ri
		const Complex J = std::polar(1.0, 2 * M_PI / 3), J2 = J * J, r1 = std::pow(r, 1.0 / 3), r12 = r1 * r1;	// r2 = J * r1, r3 = J^2 * r1;

		// Pi = P mod Ri
		Complex * const P1 = &P[0 * m / 3]; Complex * const P2 = &P[1 * m / 3]; Complex * const P3 = &P[2 * m / 3];
		for (size_t i = 0; i < m / 3; ++i)
		{
			const Complex u0 = P1[i], u1 = r1 * P2[i], u2 = r12 * P3[i];
			P1[i] = u0 +      u1 +      u2;
			P2[i] = u0 +  J * u1 + J2 * u2;
			P3[i] = u0 + J2 * u1 +  J * u2;
		}

		// Pi = P^2 mod Ri => divide and conquer
		square_fast(P1, m / 3, r1); square_fast(P2, m / 3, J * r1); square_fast(P3, m / 3, J2 * r1);

		// Chinese remainder theorem or simply the inverse of direct butterfly
		for (size_t i = 0; i < m / 3; ++i)
		{
			const Complex u0 = P1[i], u1 = P2[i], u2 = P3[i];
			P1[i] = (u0 +      u1 +      u2) /  3.0;
			P2[i] = (u0 + J2 * u1 + J  * u2) / (3.0 * r1);
			P3[i] = (u0 + J  * u1 + J2 * u2) / (3.0 * r12);
		}
	}
	else if (m % 2 == 0)
	{
		// R = R1 * R2, where Ri = x^{m/2} - ri
		const Complex r1 = std::sqrt(r);	// r2 = -r1;

		// Pi = P mod Ri: if P = a[0] + a[1].x + ... then Pi mod (x^{m/2} - ri) = (a[0] + ri * a[m/2]) + ...
		Complex * const P1 = &P[0 * m / 2]; Complex * const P2 = &P[1 * m / 2];
		for (size_t i = 0; i < m / 2; ++i)
		{
			const Complex u0 = P1[i], u1 = r1 * P2[i];
			P1[i] = u0 + u1;
			P2[i] = u0 - u1;
		}

		// P1 = P^2 mod R1, P2 = P^2 mod R2 => divide and conquer
		square_fast(P1, m / 2, r1); square_fast(P2, m / 2, -r1);

		// Chinese remainder theorem or simply the inverse of direct butterfly
		for (size_t i = 0; i < m / 2; ++i)
		{
			const Complex u0 = P1[i], u1 = P2[i];
			P1[i] = (u1 + u0) /  2.0;
			P2[i] = (u0 - u1) / (2.0 * r1);
		}
	}
}

// Q = P^2 mod x^n - r
static void square_slow(Complex * const Q, const Complex * const P, const size_t n, const Complex & r)
{
	for (size_t k = 0; k < n; ++k)
	{
		Complex l = 0; for (size_t i = 0; i <= k; ++i) l += P[i] * P[k - i];
		Complex h = 0; for (size_t i = k + 1; i < n; ++i) h += P[i] * P[k - i + n];
		Q[k] = l + r * h;
	}
}

static void display(const Complex * const P, const size_t n)
{
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
}

int main()
{
	std::srand((unsigned int)(std::time(nullptr)));

	const size_t n = 2*2*3*3;
	Complex P[n]; for (size_t i = 0; i < n; ++i) P[i] = std::rand() % 10;
	Complex r = 5;	// modulo x^n - r

	std::cout << "Mod(";
	display(P, n);
	std::cout << ", x^" << n << " - " << r.real() << ")^2" << std::endl;

	Complex Q[n]; square_slow(Q, P, n, r);
	std::cout << " = ";
	display(Q, n);
	std::cout << std::endl;

	square_fast(P, n, r);
	std::cout << " = ";
	display(P, n);
	std::cout << std::endl;

	return EXIT_SUCCESS;
}
