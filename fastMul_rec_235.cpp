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

	if (m % 2 == 0)
	{
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
	else if (m % 3 == 0)
	{
		// R = R0 * R1 * R2, where Ri = x^{m/3} - ri
		const Complex J = std::polar(1.0, 2 * M_PI / 3), J2 = J * J, r0 = std::pow(r, 1.0 / 3);	// ri = J^i * r0
		const Complex r02 = r0 * r0;

		// Pi = P mod Ri
		Complex * const P0 = &P[0 * m / 3]; Complex * const P1 = &P[1 * m / 3]; Complex * const P2 = &P[2 * m / 3];

		for (size_t i = 0; i < m / 3; ++i)
		{
			const Complex u0 = P0[i], u1 = r0 * P1[i], u2 = r02 * P2[i];
			P0[i] = u0 +      u1 +      u2;
			P1[i] = u0 + J  * u1 + J2 * u2;
			P2[i] = u0 + J2 * u1 + J  * u2;
		}

		// Pi = P^2 mod Ri => divide and conquer
		square_fast(P0, m / 3, r0); square_fast(P1, m / 3, J * r0); square_fast(P2, m / 3, J2 * r0);

		// Chinese remainder theorem or simply the inverse of direct butterfly
		for (size_t i = 0; i < m / 3; ++i)
		{
			const Complex u0 = P0[i], u1 = P1[i], u2 = P2[i];
			P0[i] = (u0 +      u1 +      u2) /  3.0;
			P1[i] = (u0 + J2 * u1 + J  * u2) / (3.0 * r0);
			P2[i] = (u0 + J  * u1 + J2 * u2) / (3.0 * r02);
		}
	}
	else if (m % 5 == 0)
	{
		// R = R0 * R1 * R2 * R3 * R4, where Ri = x^{m/3} - ri
		const Complex K = std::polar(1.0, 2 * M_PI / 5), K2 = K * K, K3 = K * K2, K4 = K2 * K2, r0 = std::pow(r, 1.0 / 5);	// ri = K^i * r0
		const Complex r02 = r0 * r0, r03 = r0 * r02, r04 = r02 * r02;

		// Pi = P mod Ri
		Complex * const P0 = &P[0 * m / 5]; Complex * const P1 = &P[1 * m / 5];
		Complex * const P2 = &P[2 * m / 5]; Complex * const P3 = &P[3 * m / 5];  Complex * const P4 = &P[4 * m / 5];

		for (size_t i = 0; i < m / 5; ++i)
		{
			const Complex u0 = P0[i], u1 = r0 * P1[i], u2 = r02 * P2[i], u3 = r03 * P3[i], u4 = r04 * P4[i];
			P0[i] = u0 +      u1 +      u2 +      u3 +      u4;
			P1[i] = u0 + K  * u1 + K2 * u2 + K3 * u3 + K4 * u4;
			P2[i] = u0 + K2 * u1 + K4 * u2 + K  * u3 + K3 * u4;
			P3[i] = u0 + K3 * u1 + K  * u2 + K4 * u3 + K2 * u4;
			P4[i] = u0 + K4 * u1 + K3 * u2 + K2 * u3 + K  * u4;
		}

		// Pi = P^2 mod Ri => divide and conquer
		square_fast(P0, m / 5, r0); square_fast(P1, m / 5, K * r0); square_fast(P2, m / 5, K2 * r0); square_fast(P3, m / 5, K3 * r0); square_fast(P4, m / 5, K4 * r0);

		// Chinese remainder theorem or simply the inverse of direct butterfly
		for (size_t i = 0; i < m / 5; ++i)
		{
			const Complex u0 = P0[i], u1 = P1[i], u2 = P2[i], u3 = P3[i], u4 = P4[i];
			P0[i] = (u0 +      u1 +      u2 +      u3 +      u4) /  5.0;
			P1[i] = (u0 + K4 * u1 + K3 * u2 + K2 * u3 + K  * u4) / (5.0 * r0);
			P2[i] = (u0 + K3 * u1 + K  * u2 + K4 * u3 + K2 * u4) / (5.0 * r02);
			P3[i] = (u0 + K2 * u1 + K4 * u2 + K  * u3 + K3 * u4) / (5.0 * r03);
			P4[i] = (u0 + K  * u1 + K2 * u2 + K3 * u3 + K4 * u4) / (5.0 * r04);
		}
	}
}

// Q = P^2 mod x^n - r
static void square_slow(Complex * const Q, const Complex * const P, const size_t n, const Complex & r)
{
	for (size_t j = 0; j < n; ++j)
	{
		Complex l = 0; for (size_t i = 0; i <= j; ++i) l += P[i] * P[j - i];
		Complex h = 0; for (size_t i = j + 1; i < n; ++i) h += P[i] * P[j - i + n];
		Q[j] = l + r * h;
	}
}

inline void check(const Complex * const P, const Complex * const Q, const size_t n)
{
	bool error = false;
	long long a_max = 0;
	for (size_t i = 0; i < n; ++i)
	{
		const long long a = std::llrint(P[i].real()), b = std::llrint(Q[i].real());
		if (a != b) { error = true; std::cout << i << ": " << a << " != " << b << std::endl; }
		a_max = std::max(a_max, a);
	}
	std::cout << (error ? "N" : "") << "OK, max = " << a_max << std::endl;
}

inline void display(const Complex * const P, const size_t n)
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

	const size_t n = 2*2*3*3*5*5;
	Complex P[n]; for (size_t i = 0; i < n; ++i) P[i] = std::rand() % 1000;
	Complex r = 7;	// modulo x^n - r

	// std::cout << "Mod(";
	// display(P, n);
	// std::cout << ", x^" << n << " - " << r.real() << ")^2" << std::endl;

	Complex Q[n]; square_slow(Q, P, n, r);
	// std::cout << " = ";
	// display(Q, n);
	// std::cout << std::endl;

	square_fast(P, n, r);
	// std::cout << " = ";
	// display(P, n);
	// std::cout << std::endl;

	check(P, Q, n);

	return EXIT_SUCCESS;
}
