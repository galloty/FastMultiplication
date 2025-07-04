/*
Copyright 2024, Yves Gallot

FastMultiplication is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <iostream>
#include <cstdlib>
#include <ctime>

// a finite field
class GF
{
private:
	uint32_t n;

public:
	// Let p be the order of the finite field. The convolution is modulo x^m - r then we must have:
	// - p > m * (b - 1)^2, where m is the transform size and b the base representation.
	// - a mth root of r must exist.
	// If b = 1000, m = 2^2 * 3^2 * 5^2 and r = 7, the first condition is p > 898200900.
	// If p = 913262401 then 243771734^m = 7 (mod p).

	static const uint32_t _p = 913262401u;
	static const uint32_t _primroot = 11u;
	static const uint32_t _root1 = 422058578u;	// 11 is a primitive root. 11^{(p - 1)/m} = 422058578 then 422058578^m = 1
	static const uint32_t _root7 = 243771734u;

public:
	GF() {}
	GF(const uint32_t l) : n(l) {}

	uint32_t get() const { return n; }

	GF & operator+=(const GF & rhs) { const uint32_t c = (n >= _p - rhs.n) ? _p : 0; n += rhs.n; n -= c; return *this; }
	GF & operator-=(const GF & rhs) { const uint32_t c = (n < rhs.n) ? _p : 0; n -= rhs.n; n += c; return *this; }
	GF & operator*=(const GF & rhs) { *this = uint32_t((n * uint64_t(rhs.n)) % _p); return *this; }

	GF operator+(const GF & rhs) const { GF r = *this; r += rhs; return r; }
	GF operator-(const GF & rhs) const { GF r = *this; r -= rhs; return r; }
	GF operator*(const GF & rhs) const { GF r = *this; r *= rhs; return r; }

	GF pow(const uint32_t e) const
	{
		if (e == 0) return GF(1u);

		GF r = GF(1u), y = *this;
		for (uint32_t i = e; i != 1; i /= 2)
		{
			if (i % 2 != 0) r *= y;
			y *= y;
		}
		r *= y;

		return r;
	}

	GF invert() const { return pow(_p - 2); }

	static const GF root_1(const uint32_t n) { return GF(_root1).pow(n); }
	static const GF root_7(const uint32_t n) { return GF(_root7).pow(n); }
	static const GF root_nth(const uint32_t n) { return GF(_primroot).pow((_p - 1) / n); }
};

// P = P^2 mod R, where R = x^m - r
static void square_fast(GF * const P, const size_t j, const size_t m, const size_t s)
{
	if (m == 1) { P[0] *= P[0]; return; } 	// P is a scalar

	if (m % 2 == 0)
	{
		// R = R0 * R1, where Ri = x^{m/2} - ri
		const GF r0 = GF::root_7(m / 2) * GF::root_1(j * m / 2);	// r1 = -r0;
 
		// Pi = P mod Ri: if P = a[0] + a[1].x + ... then Pi mod (x^{m/2} - ri) = (a[0] + ri * a[m/2]) + ...
		GF * const P0 = &P[0 * m / 2]; GF * const P1 = &P[1 * m / 2];

		for (size_t i = 0; i < m / 2; ++i)
		{
			const GF u0 = P0[i], u1 = r0 * P1[i];
			P0[i] = u0 + u1;
			P1[i] = u0 - u1;
		}

		// R1 = x^{m/2} + r0, P1 must be twisted such that R1' = x^{m/2} - r0
		for (size_t i = 0; i < m / 2; ++i) P1[i] *= GF::root_1(s * i);

		// P1 = P^2 mod R1, P2 = P^2 mod R2 => divide and conquer
		square_fast(P0, j + 0 * s, m / 2, s * 2); square_fast(P1, j /*+ 1 * s*/, m / 2, s * 2);

		for (size_t i = 0; i < m / 2; ++i) P1[i] *= GF::root_1(s * i).invert();	// untwist

		// Chinese remainder theorem or simply the inverse of direct butterfly
		const GF inv2 = GF(2u).invert(), invr0 = r0.invert();
		for (size_t i = 0; i < m / 2; ++i)
		{
			const GF u0 = P0[i], u1 = P1[i];
			P0[i] = (u1 + u0) * inv2;
			P1[i] = (u0 - u1) * inv2 * invr0;
		}
	}
	else if (m % 3 == 0)
	{
		// R = R0 * R1 * R2, where Ri = x^{m/3} - ri
		static const GF J = GF::root_nth(3u), J2 = J * J;
		const GF r0 = GF::root_7(m / 3) * GF::root_1(j * m / 3);	// ri = J^i * r0
		const GF r02 = r0 * r0;

		// Pi = P mod Ri
		GF * const P0 = &P[0 * m / 3]; GF * const P1 = &P[1 * m / 3]; GF * const P2 = &P[2 * m / 3];

		for (size_t i = 0; i < m / 3; ++i)
		{
			const GF u0 = P0[i], u1 = r0 * P1[i], u2 = r02 * P2[i];
			P0[i] = u0 +      u1 +      u2;
			P1[i] = u0 + J  * u1 + J2 * u2;
			P2[i] = u0 + J2 * u1 + J  * u2;
		}

		// P1 and P2 are twisted such that R1' = R2' = x^{m/3} - r0
		for (size_t i = 0; i < m / 3; ++i) P1[i] *= GF::root_1(1 * s * i);
		for (size_t i = 0; i < m / 3; ++i) P2[i] *= GF::root_1(2 * s * i);

		// Pi = P^2 mod Ri => divide and conquer
		square_fast(P0, j + 0 * s, m / 3, s * 3); square_fast(P1, j /*+ 1 * s*/, m / 3, s * 3); square_fast(P2, j /*+ 2 * s*/, m / 3, s * 3);

		for (size_t i = 0; i < m / 3; ++i) P1[i] *= GF::root_1(1 * s * i).invert();	// untwist
		for (size_t i = 0; i < m / 3; ++i) P2[i] *= GF::root_1(2 * s * i).invert();	;

		// Chinese remainder theorem or simply the inverse of direct butterfly
		const GF inv3 = GF(3u).invert(), invr0 = r0.invert(), invr02 = invr0 * invr0;
		for (size_t i = 0; i < m / 3; ++i)
		{
			const GF u0 = P0[i], u1 = P1[i], u2 = P2[i];
			P0[i] = (u0 +      u1 +      u2) * inv3;
			P1[i] = (u0 + J2 * u1 + J  * u2) * inv3 * invr0;
			P2[i] = (u0 + J  * u1 + J2 * u2) * inv3 * invr02;
		}
	}
	else if (m % 5 == 0)
	{
		// R = R0 * R1 * R2 * R3 * R4, where Ri = x^{m/3} - ri
		static const GF K = GF::root_nth(5u), K2 = K * K, K3 = K * K2, K4 = K2 * K2;
		const GF r0 = GF::root_7(m / 5) * GF::root_1(j * m / 5);	// ri = K^i * r0
		const GF r02 = r0 * r0, r03 = r0 * r02, r04 = r02 * r02;

		// Pi = P mod Ri
		GF * const P0 = &P[0 * m / 5]; GF * const P1 = &P[1 * m / 5];
		GF * const P2 = &P[2 * m / 5]; GF * const P3 = &P[3 * m / 5];  GF * const P4 = &P[4 * m / 5];

		for (size_t i = 0; i < m / 5; ++i)
		{
			const GF u0 = P0[i], u1 = r0 * P1[i], u2 = r02 * P2[i], u3 = r03 * P3[i], u4 = r04 * P4[i];
			P0[i] = u0 +      u1 +      u2 +      u3 +      u4;
			P1[i] = u0 + K  * u1 + K2 * u2 + K3 * u3 + K4 * u4;
			P2[i] = u0 + K2 * u1 + K4 * u2 + K  * u3 + K3 * u4;
			P3[i] = u0 + K3 * u1 + K  * u2 + K4 * u3 + K2 * u4;
			P4[i] = u0 + K4 * u1 + K3 * u2 + K2 * u3 + K  * u4;
		}

		// Pi = P^2 mod Ri => divide and conquer
		square_fast(P0, j + 0 * s, m / 5, s * 5); square_fast(P1, j + 1 * s, m / 5, s * 5); square_fast(P2, j + 2 * s, m / 5, s * 5);
		square_fast(P3, j + 3 * s, m / 5, s * 5); square_fast(P4, j + 4 * s, m / 5, s * 5);

		// Chinese remainder theorem or simply the inverse of direct butterfly
		const GF inv5 = GF(5u).invert(), invr0 = r0.invert(), invr02 = invr0 * invr0, invr03 = invr0 * invr02, invr04 = invr02 * invr02;
		for (size_t i = 0; i < m / 5; ++i)
		{
			const GF u0 = P0[i], u1 = P1[i], u2 = P2[i], u3 = P3[i], u4 = P4[i];
			P0[i] = (u0 +      u1 +      u2 +      u3 +      u4) * inv5;
			P1[i] = (u0 + K4 * u1 + K3 * u2 + K2 * u3 + K  * u4) * inv5 * invr0;
			P2[i] = (u0 + K3 * u1 + K  * u2 + K4 * u3 + K2 * u4) * inv5 * invr02;
			P3[i] = (u0 + K2 * u1 + K4 * u2 + K  * u3 + K3 * u4) * inv5 * invr03;
			P4[i] = (u0 + K  * u1 + K2 * u2 + K3 * u3 + K4 * u4) * inv5 * invr04;
		}
	}
}

// Q = P^2 mod x^n - r
static void square_slow(GF * const Q, const GF * const P, const size_t n, const GF & r)
{
	for (size_t j = 0; j < n; ++j)
	{
		GF l = 0; for (size_t i = 0; i <= j; ++i) l += P[i] * P[j - i];
		GF h = 0; for (size_t i = j + 1; i < n; ++i) h += P[i] * P[j - i + n];
		Q[j] = l + r * h;
	}
}

inline void check(const GF * const P, const GF * const Q, const size_t n)
{
	bool error = false;
	uint32_t a_max = 0;
	for (size_t i = 0; i < n; ++i)
	{
		const uint32_t a = P[i].get(), b = Q[i].get();
		if (a != b) { error = true; std::cout << i << ": " << a << " != " << b << std::endl; }
		a_max = std::max(a_max, a);
	}
	std::cout << (error ? "N" : "") << "OK, max = " << a_max << std::endl;
}

inline void display(const GF * const P, const size_t n)
{
	for (size_t i = n - 1, j = 0; j < n; --i, ++j)
	{
		std::cout << P[i].get();
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
	GF P[n]; for (size_t i = 0; i < n; ++i) P[i] = GF(uint32_t(std::rand()) % 1000u);
	const GF r = GF(7u);	// modulo x^n - r

	// std::cout << "Mod(";
	// display(P, n);
	// std::cout << ", x^" << n << " - " << r.get() << ")^2" << std::endl;

	GF Q[n]; square_slow(Q, P, n, r);
	// std::cout << " = ";
	// display(Q, n);
	// std::cout << std::endl;

	square_fast(P, 0, n, 1);	// s = n / m
	// std::cout << " = ";
	// display(P, n);
	// std::cout << std::endl;

	check(P, Q, n);

	return EXIT_SUCCESS;
}
