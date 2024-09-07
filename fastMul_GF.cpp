/*
Copyright 2024, Yves Gallot

FastMultiplication is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <iostream>

// Modulo 2^61 - 1
class Zp
{
private:
	static const uint64_t _p = (uint64_t(1) << 61) - 1;
	uint64_t _n;

	static uint64_t add(const uint64_t a, const uint64_t b) { return a + b - ((a >= _p - b) ? _p : 0); }
	static uint64_t sub(const uint64_t a, const uint64_t b) { return a - b + ((a < b) ? _p : 0); }
	static uint64_t mul(const uint64_t a, const uint64_t b) { return uint64_t((a * __uint128_t(b)) % _p); }

public:
	Zp() {}
	Zp(const uint64_t n) : _n(n) {}

	uint64_t get() const { return _n; }
	int64_t get_int() const { return (_n >= _p / 2) ? int64_t(_n - _p) : int64_t(_n); }

	Zp operator-() const { return Zp((_n == 0) ? 0 : _p - _n); }

	Zp operator+(const Zp & rhs) const { return Zp(add(_n, rhs._n)); }
	Zp operator-(const Zp & rhs) const { return Zp(sub(_n, rhs._n)); }
	Zp operator*(const Zp & rhs) const { return Zp(mul(_n, rhs._n)); }

	Zp & operator*=(const Zp & rhs) { _n = mul(_n, rhs._n); return *this; }

	Zp pow(const uint64_t e) const
	{
		if (e == 0) return Zp(1);
		Zp r = Zp(1), y = *this;
		for (uint64_t i = e; i != 1; i /= 2) { if (i % 2 != 0) r *= y; y *= y; }
		return r * y;
	}

	Zp invert() const { return pow(_p - 2); }
};

// GF((2^61 - 1)^2)
class GF
{
private:
	Zp _a, _b;
	static const uint64_t _h_a = 2147483648ull, _h_b = 1272521237944691271ull;	// a primitive root of order 2^62 which is a root of (0, 1).

	GF(const Zp & a, const Zp & b) : _a(a), _b(b) {}

public:
	GF() {}
	GF(const int64_t a, const int64_t b) : _a(uint64_t(a)), _b(uint64_t(b)) {}	// a, b >= 0

	const Zp & a() const { return _a; }
	const Zp & b() const { return _b; }

	GF conj() const { return GF(_a, -_b); }

	GF operator+(const GF & rhs) const { return GF(_a + rhs._a, _b + rhs._b); }
	GF operator-(const GF & rhs) const { return GF(_a - rhs._a, _b - rhs._b); }
	GF operator*(const GF & rhs) const { return GF(_a * rhs._a - _b * rhs._b, _a * rhs._b + _b * rhs._a); }

	GF & operator*=(const GF & rhs) { *this = *this * rhs; return *this; }

	// GF invert() const
	// {
	// 	const Zp t = Zp(_a * _a + _b * _b).invert();
	// 	return GF(_a * t, (-_b) * t);
	// }

	GF pow(const uint64_t e) const
	{
		if (e == 0) return GF(1, 0);
		GF r = GF(1, 0), y = *this;
		for (uint64_t i = e; i != 1; i /= 2) { if (i % 2 != 0) r *= y; y *= y; }
		return r * y;
	}

	static const GF primroot_n(const uint64_t n) { return GF(Zp(_h_a), Zp(_h_b)).pow((uint64_t(1) << 62) / n); }
};

constexpr size_t bitrev(const size_t i, const size_t n)
{
	size_t r = 0; for (size_t k = n, j = i; k > 1; k /= 2, j /= 2) r = (2 * r) | (j % 2); return r;
}

static void roots(GF * const w, const size_t n)
{
	// Note that:
	// - w[1] = 2^30 * (1, 1)
	// - if w[2k] = (a, b) then w[2k + 1] = (-b, a)
	for (size_t s = 1; s < n; s *= 2)
	{
		const GF r_s = GF::primroot_n(8 * s);
		for (size_t j = 0; j < s; ++j)
		{
			w[s + j] = r_s.pow(bitrev(j, 4 * s) + 1);
			// std::cout << s << ", " << j << ": " << w[s + j].a().get_int() << ", " << w[s + j].b().get_int() << std::endl;
		}
	}
}

static void forward2(GF * const x, const GF * const w, const size_t n)
{
	for (size_t s = 1, m = n / 2; m >= 1; s *= 2, m /= 2)
	{
		for (size_t j = 0; j < s; ++j)
		{
			const GF w_sj = w[s + j];

			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = j * 2 * m + i;
				const GF u0 = x[k + 0 * m], u1 = x[k + 1 * m] * w_sj;
				x[k + 0 * m] = u0 + u1; x[k + 1 * m] = u0 - u1;
			}
		}
	}
}

static void backward2(GF * const x, const GF * const w, const size_t n)
{
	for (size_t s = n / 2, m = 1; s >= 1; s /= 2, m *= 2)
	{
		for (size_t j = 0; j < s; ++j)
		{
			const GF w_sj_inv = w[s + j].conj();	// conjugate is inverse

			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = j * 2 * m + i;
				const GF u0 = x[k + 0 * m], u1 = x[k + 1 * m];
				x[k + 0 * m] = u0 + u1; x[k + 1 * m] = (u0 - u1) * w_sj_inv;
			}
		}
	}
}

static void square(GF * const x, const size_t n)
{
	for (size_t k = 0; k < n; ++k) x[k] *= x[k];
}

int main()
{
	const size_t n = 16;

	// (25M * (n - 1) + 1)^2 * n = 2.25e18 < 2^61 - 1 = 2.30e18
	int64_t a[n]; for (size_t i = 0; i < n; ++i) a[i] = 25000000*i + 1;

	GF w[n / 2]; roots(w, n / 2);

	// A transform of length n / 2 can be used because the polynomial of GF is X^2 + 1
	GF x[n / 2]; for (size_t i = 0; i < n / 2; ++i) x[i] = GF(a[i], a[i + n / 2]);

	// Can be optimized using Radix-4 which has 3/4 as many GF multiplies as Radix-2
	forward2(x, w, n / 2);
	square(x, n / 2);
	backward2(x, w, n / 2);

	std::cout << "Mod(1 + 25000001*x + 50000001*x^2 + 75000001*x^3 + 100000001*x^4 + 125000001*x^5 + 150000001*x^6 + 175000001*x^7 + 200000001*x^8 + 225000001*x^9 + 250000001*x^10 + 275000001*x^11 + 300000001*x^12 + 325000001*x^13 + 350000001*x^14 + 375000001*x^15, x^16 + 1)^2" << std::endl;
	std::cout << " = 350000006000000016*x^15 + 143750004500000014*x^14 - 34999996899999988*x^13 - 187499998199999990*x^12 - 314999999399999992*x^11 - 418750000499999994*x^10 - 500000001499999996*x^9 - 560000002399999998*x^8 - 600000003200000000*x^7 - 621250003900000002*x^6 - 625000004500000004*x^5 - 612500005000000006*x^4 - 585000005400000008*x^3 - 543750005700000010*x^2 - 490000005900000012*x - 425000006000000014" << std::endl << "   ";

	const Zp r = Zp(n / 2).invert();
	for (size_t i = n - 1, j = 0; j < n; --i, ++j)
	{
		const auto a_j = (i >= n / 2) ? x[i - n / 2].b() : x[i].a();
		const int64_t k = Zp(a_j * r).get_int();

		if (j != 0) std::cout << ((k >= 0) ? " + ": " - ");
		std::cout << std::abs(k);
		if (i > 1) std::cout << "*x^" << i; else if (i > 0) std::cout << "*x";
	} 
	std::cout << std::endl;

	return EXIT_SUCCESS;
}
