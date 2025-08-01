/*
Copyright 2024, Yves Gallot

FastMultiplication is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <iostream>
#include <cstdlib>
#include <ctime>

// A prime finite field such that p = k * n + 1 and a nth root of r exists
class GF
{
private:
	uint32_t _n;

public:
	// Let p be the order of the finite field. We must have:
	// - p > n * (b - 1)^2, where n is the transform size and b the base representation.
	// - a nth root of r must exist.
	// If b = 1000, n = 2^2 * 3^2 * 5^2 and r = 7, the first condition is p > 898200900.
	// If p = 913262401 then 243771734^m = 7 (mod p).

	static const uint32_t _p = 913262401u;
	static const uint32_t _primroot = 11u;
	static const uint32_t _root1 = 422058578u;	// 422058578^n = 1
	static const uint32_t _root7 = 243771734u;	// 243771734^n = 7

public:
	GF() {}
	GF(const uint32_t n) : _n(n) {}

	uint32_t get() const { return _n; }

	GF & operator+=(const GF & rhs) { const uint32_t c = (_n >= _p - rhs._n) ? _p : 0; _n += rhs._n; _n -= c; return *this; }
	GF & operator-=(const GF & rhs) { const uint32_t c = (_n < rhs._n) ? _p : 0; _n -= rhs._n; _n += c; return *this; }
	GF & operator*=(const GF & rhs) { *this = uint32_t((_n * uint64_t(rhs._n)) % _p); return *this; }

	GF operator+(const GF & rhs) const { GF r = *this; r += rhs; return r; }
	GF operator-(const GF & rhs) const { GF r = *this; r -= rhs; return r; }
	GF operator*(const GF & rhs) const { GF r = *this; r *= rhs; return r; }

	GF half() const { return GF((_n % 2 == 0) ? _n / 2 : ((_n - 1) / 2 + (_p + 1) / 2)); }

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

	static const GF root_900th_1() { return GF(_root1); }
	static const GF root_900th_7() { return GF(_root7); }
	static GF root_nth(const uint32_t n) { return GF(_primroot).pow((_p - 1) / n); }
};

// P = P^2 mod x^n - r
class conv_fast
{
private:
	const size_t _n;
	GF * const _root;
	GF * const _invroot;
	const GF J = GF::root_nth(3u);	// Radix-3
	const GF K = GF::root_nth(5u), K1 = K + GF(1u), K2 = K * K;	// Radix-5
	const GF F2 = (K1 * K2).half(), F3 = K * K1, F4 = K * (K2 + GF(1u)), F5 = (F3 + F4).half();

	// Generalization of bit-reversal permutation
	static constexpr size_t reversal(const size_t i, const size_t n)
	{
		size_t r = 0, k = n, j = i;
		// The order should be the inverse of the order of the square function
		for (; k % 2 == 0; k /= 2, j /= 2) r = 2 * r + j % 2;
		for (; k % 3 == 0; k /= 3, j /= 3) r = 3 * r + j % 3;
		for (; k % 5 == 0; k /= 5, j /= 5) r = 5 * r + j % 5;
		return r;
	}

	void set_roots(const size_t s, const size_t m) const
	{
		GF * const root = _root;
		GF * const invroot = _invroot;
		for (size_t j = 0; j < s; ++j)
		{
			const GF r = (GF::root_900th_7() * GF::root_900th_1().pow(reversal(j, s))).pow(m);
			root[s + j] = r; invroot[s + j] = r.invert();
		}
	}

public:
	conv_fast(const size_t n) : _n(n), _root(new GF[n]), _invroot(new GF[n])
	{
		// Init precomputed tables
		size_t m = _n, s = 1;
		while (m != 1)
		{
			if      (m % 5 == 0) { m /= 5; set_roots(s, m); s *= 5; }
			else if (m % 3 == 0) { m /= 3; set_roots(s, m); s *= 3; }
			else if (m % 2 == 0) { m /= 2; set_roots(s, m); s *= 2; }
		}
	}

	~conv_fast()
	{
		delete[] _root;
		delete[] _invroot;
	}

	void square(GF * const P) const
	{
		const size_t n = _n;
		const GF * const root = _root;
		const GF * const invroot = _invroot;

		size_t m = n, s = 1;
		while (m != 1)
		{
			if (m % 5 == 0)
			{
				m /= 5;

				for (size_t j = 0; j < s; ++j)
				{
					const GF r = root[s + j], r2 = r * r, r3 = r * r2, r4 = r2 * r2;

					for (size_t i = 0; i < m; ++i)
					{
						const size_t k = 5 * m * j + i;
						const GF u0 = P[k + 0 * m], u1 = r * P[k + 1 * m], u2 = r2 * P[k + 2 * m], u3 = r3 * P[k + 3 * m], u4 = r4 * P[k + 4 * m];

						// const GF K3 = K * K2, K4 = K2 * K2;
						// P[k + 0 * m] = u0 +      u1 +      u2 +      u3 +      u4;
						// P[k + 1 * m] = u0 + K  * u1 + K2 * u2 + K3 * u3 + K4 * u4;
						// P[k + 2 * m] = u0 + K2 * u1 + K4 * u2 + K  * u3 + K3 * u4;
						// P[k + 3 * m] = u0 + K3 * u1 + K  * u2 + K4 * u3 + K2 * u4;
						// P[k + 4 * m] = u0 + K4 * u1 + K3 * u2 + K2 * u3 + K  * u4;

						const GF v1 = u1 + u4, v4 = u1 - u4, v2 = u2 + u3, v3 = u2 - u3;
						const GF z0 = v1 + v2;
						const GF x1 = u0, x2 = (v1 - v2) * F2, x3 = v3 * F3, x4 = v4 * F4, x5 = (v4 - v3) * F5;
						const GF y1 = x1 + x2, y2 = x1 - x2, y3 = x3 + x5, y4 = x4 - x5;
						const GF z1 = y1 + y4, z4 = y1 - y4, z2 = y2 + y3, z3 = y2 - y3;
						P[k + 0 * m] = z0 + u0; P[k + 1 * m] = z2 - u4; P[k + 2 * m] = z4 - u2; P[k + 3 * m] = z1 - u3; P[k + 4 * m] = z3 - u1;
					}
				}

				s *= 5;
			}
			else if (m % 3 == 0)
			{
				m /= 3;

				for (size_t j = 0; j < s; ++j)
				{
					const GF r = root[s + j], r2 = r * r;

					for (size_t i = 0; i < m; ++i)
					{
						const size_t k = 3 * m * j + i;
						const GF u0 = P[k + 0 * m], u1 = r * P[k + 1 * m], u2 = r2 * P[k + 2 * m];
						const GF t = J * (u1 - u2);	// J^2 = -(1 + J)
						P[k + 0 * m] = u0 + u1 + u2;
						P[k + 1 * m] = u0 - u2 + t;
						P[k + 2 * m] = u0 - u1 - t;
					}
				}

				s *= 3;
			}
			else if (m % 2 == 0)
			{
				m /= 2;

				for (size_t j = 0; j < s; ++j)
				{
					const GF r = root[s + j];

					for (size_t i = 0; i < m; ++i)
					{
						const size_t k = 2 * m * j + i;
						const GF u0 = P[k + 0 * m], u1 = r * P[k + 1 * m];
						P[k + 0 * m] = u0 + u1;
						P[k + 1 * m] = u0 - u1;
					}
				}

				s *= 2;
			}
		}

		for (size_t k = 0; k < n; ++k) P[k] *= P[k];

		while (s != 1)
		{
			if (s % 2 == 0)
			{
				s /= 2;

				for (size_t j = 0; j < s; ++j)
				{
					const GF invr = invroot[s + j];

					for (size_t i = 0; i < m; ++i)
					{
						const size_t k = 2 * m * j + i;
						const GF u0 = P[k + 0 * m], u1 = P[k + 1 * m];
						P[k + 0 * m] =  u1 + u0;
						P[k + 1 * m] = (u0 - u1) * invr;
					}
				}

				m *= 2;
			}
			else if (s % 3 == 0)
			{
				s /= 3;

				for (size_t j = 0; j < s; ++j)
				{
					const GF invr = invroot[s + j], invr2 = invr * invr;

					for (size_t i = 0; i < m; ++i)
					{
						const size_t k = 3 * m * j + i;
						const GF u0 = P[k + 0 * m], u1 = P[k + 1 * m], u2 = P[k + 2 * m];
						const GF t = J * (u1 - u2);	// J^2 = -(1 + J)
						P[k + 0 * m] =  u0 + u1 + u2;
						P[k + 1 * m] = (u0 - u1 - t) * invr;
						P[k + 2 * m] = (u0 - u2 + t) * invr2;
					}
				}

				m *= 3;
			}
			else if (s % 5 == 0)
			{
				s /= 5;

				for (size_t j = 0; j < s; ++j)
				{
					const GF invr = invroot[s + j], invr2 = invr * invr, invr3 = invr * invr2, invr4 = invr2 * invr2;

					for (size_t i = 0; i < m; ++i)
					{
						const size_t k = 5 * m * j + i;
						const GF u0 = P[k + 0 * m], u4 = P[k + 1 * m], u3 = P[k + 2 * m], u2 = P[k + 3 * m], u1 = P[k + 4 * m];

						// const GF K3 = K * K2, K4 = K2 * K2;
						// P[k + 0 * m] =  u0 +      u1 +      u2 +      u3 +      u4;
						// P[k + 1 * m] = (u0 + K  * u1 + K2 * u2 + K3 * u3 + K4 * u4) * invr;
						// P[k + 2 * m] = (u0 + K2 * u1 + K4 * u2 + K  * u3 + K3 * u4) * invr2;
						// P[k + 3 * m] = (u0 + K3 * u1 + K  * u2 + K4 * u3 + K2 * u4) * invr3;
						// P[k + 4 * m] = (u0 + K4 * u1 + K3 * u2 + K2 * u3 + K  * u4) * invr4;

						const GF v1 = u1 + u4, v4 = u1 - u4, v2 = u2 + u3, v3 = u2 - u3;
						const GF z0 = v1 + v2;
						const GF x1 = u0, x2 = (v1 - v2) * F2, x3 = v3 * F3, x4 = v4 * F4, x5 = (v4 - v3) * F5;
						const GF y1 = x1 + x2, y2 = x1 - x2, y3 = x3 + x5, y4 = x4 - x5;
						const GF z1 = y1 + y4, z4 = y1 - y4, z2 = y2 + y3, z3 = y2 - y3;
						P[k + 0 * m] =  z0 + u0;
						P[k + 1 * m] = (z2 - u4) * invr;;
						P[k + 2 * m] = (z4 - u2) * invr2;
						P[k + 3 * m] = (z1 - u3) * invr3;
						P[k + 4 * m] = (z3 - u1) * invr4;
					}
				}

				m *= 5;
			}
		}

		const GF invn = GF(uint32_t(n)).invert();
		for (size_t k = 0; k < n; ++k) P[k] *= invn;
	}
};

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

	conv_fast conv(n);
	conv.square(P);
	// std::cout << " = ";
	// display(P, n);
	// std::cout << std::endl;

	check(P, Q, n);

	return EXIT_SUCCESS;
}
