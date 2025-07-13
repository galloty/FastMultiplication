/*
Copyright 2024, Yves Gallot

FastMultiplication is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <vector>
#include <iostream>
#include <cmath>

#include <gmp.h>

// Z/{2^61 - 1}Z: the prime field of order p = 2^61 - 1
class Z61
{
private:
	static const uint64_t _p = (uint64_t(1) << 61) - 1;
	uint64_t _n;	// 0 <= n < p

	static uint64_t add(const uint64_t a, const uint64_t b) { return a + b - ((a >= _p - b) ? _p : 0); }
	static uint64_t sub(const uint64_t a, const uint64_t b) { return a - b + ((a < b) ? _p : 0); }
	static uint64_t mul(const uint64_t a, const uint64_t b) { const __uint128_t t = a * __uint128_t(b); return add(uint64_t(t) & _p, uint64_t(t >> 61)); }
	static uint64_t lshift(const uint64_t a, const int s) { const __uint128_t t = __uint128_t(a) << s; return add(uint64_t(t) & _p, uint64_t(t >> 61)); }

public:
	Z61() {}
	explicit Z61(const uint64_t n) : _n(n) {}

	static uint64_t get_p() { return _p; }

	uint64_t get() const { return _n; }

	int64_t get_int() const { return (_n >= _p / 2) ? int64_t(_n - _p) : int64_t(_n); }
	Z61 & set_int(const int64_t i) { _n = (i < 0) ? uint64_t(i) + _p : uint64_t(i); return *this; }

	Z61 operator+(const Z61 & rhs) const { return Z61(add(_n, rhs._n)); }
	Z61 operator-(const Z61 & rhs) const { return Z61(sub(_n, rhs._n)); }
	Z61 operator*(const Z61 & rhs) const { return Z61(mul(_n, rhs._n)); }
	Z61 operator<<(const int s) const { return Z61(lshift(_n, s)); }
};

// GF((2^61 - 1)^2): the field of order p^2, p = 2^61 - 1
class GF61
{
private:
	// The irreducible polynomial is X^2 + 1: we have (0, 1)^2 = (-1, 0)
	Z61 _n0, _n1;
	// a primitive root of order 2^62 which is a root of (0, 1).
	static const uint64_t _h_order = uint64_t(1) << 62;
	static const uint64_t _h_0 = 264036120304204ull, _h_1 = 4677669021635377ull;

	explicit GF61(const Z61 & n0, const Z61 & n1) : _n0(n0), _n1(n1) {}

public:
	GF61() {}
	explicit GF61(const uint64_t n) : _n0(n), _n1(0) {}

	void get_int(int64_t & i0, int64_t & i1) const { i0 = _n0.get_int(); i1 = _n1.get_int(); }
	GF61 & set_int(const int64_t i0, const int64_t i1) { _n0.set_int(i0); _n1.set_int(i1); return *this; }

	GF61 operator+(const GF61 & rhs) const { return GF61(_n0 + rhs._n0, _n1 + rhs._n1); }
	GF61 operator-(const GF61 & rhs) const { return GF61(_n0 - rhs._n0, _n1 - rhs._n1); }
	GF61 addi(const GF61 & rhs) const { return GF61(_n0 - rhs._n1, _n1 + rhs._n0); }
	GF61 subi(const GF61 & rhs) const { return GF61(_n0 + rhs._n1, _n1 - rhs._n0); }

	GF61 sqr() const { const Z61 t = _n0 * _n1; return GF61(_n0 * _n0 -  _n1 * _n1, t + t); }
	GF61 mul(const GF61 & rhs) const { return GF61(_n0 * rhs._n0 -  _n1 * rhs._n1, _n1 * rhs._n0 + _n0 * rhs._n1); }
	GF61 mulconj(const GF61 & rhs) const { return GF61(_n0 * rhs._n0 + _n1 * rhs._n1, _n1 * rhs._n0 - _n0 * rhs._n1); }

	GF61 operator<<(const int s) const { return GF61(_n0 << s, _n1 << s); }

	GF61 pow(const uint64_t e) const
	{
		if (e == 0) return GF61(1u);
		GF61 r = GF61(1u), y = *this;
		for (uint64_t i = e; i != 1; i /= 2) { if (i % 2 != 0) r = r.mul(y); y = y.sqr(); }
		return r.mul(y);
	}

	static const GF61 primroot_n(const size_t n) { return GF61(Z61(_h_0), Z61(_h_1)).pow(_h_order / n); }
};

class Power
{
private:
	std::vector<GF61> _vz;
	std::vector<GF61> _vr;
	const int _snorm61;
	int32_t _base;

private:
	// Bit-reversal permutation
	static constexpr size_t bitrev(const size_t i, const size_t n)
	{
		size_t r = 0;
		for (size_t k = n, j = i; k != 1; k /= 2, j /= 2) r = (2 * r) | (j % 2);
		return r;
	}

	// Radix-4: multiplication by (0, 1) is free then 3 multiplications are needed
	void forward4(const size_t m, const size_t s)
	{
		GF61 * const z = _vz.data();
		const GF61 * const r = _vr.data();

		for (size_t j = 0; j < s; ++j)
		{
			const GF61 r1 = r[s + j], r2 = r[2 * (s + j)], r3 = r1.mul(r2);

			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = 4 * m * j + i;
				const GF61 u0 = z[k + 0 * m], u1 = z[k + 1 * m].mul(r2), u2 = z[k + 2 * m].mul(r1), u3 = z[k + 3 * m].mul(r3);
				const GF61 v0 = u0 + u2, v1 = u1 + u3, v2 = u0 - u2, v3 = u1 - u3;
				z[k + 0 * m] = v0 + v1; z[k + 1 * m] = v0 - v1;
				z[k + 2 * m] = v2.addi(v3); z[k + 3 * m] = v2.subi(v3);
			}
		}
	}

	// Inverse Radix-4: 3 multiplications
	void backward4(const size_t m, const size_t s)
	{
		GF61 * const z = _vz.data();
		const GF61 * const r = _vr.data();

		for (size_t j = 0; j < s; ++j)
		{
			const GF61 r1 = r[s + j], r2 = r[2 * (s + j)], r3 = r1.mul(r2);

			for (size_t i = 0; i < m; ++i)
			{
				const size_t k = 4 * m * j + i;
				const GF61 u0 = z[k + 0 * m], u1 = z[k + 1 * m], u2 = z[k + 2 * m], u3 = z[k + 3 * m];
				const GF61 v0 = u0 + u1, v1 = u0 - u1, v2 = u2 + u3, v3 = u3 - u2;
				z[k + 0 * m] = v0 + v2; z[k + 2 * m] = (v0 - v2).mulconj(r1);
				z[k + 1 * m] = v1.addi(v3).mulconj(r2); z[k + 3 * m] = v1.subi(v3).mulconj(r3);
			}
		}
	}

	// Radix-2, pointwise square, inverse radix-2
	void square2x2()
	{
		const size_t n = _vz.size();
		GF61 * const z = _vz.data();
		const GF61 * const r = _vr.data();

		for (size_t j = 0; j < n / 4; ++j)
		{
			const GF61 r0 = r[n / 4 + j];

			const size_t k = 4 * j;
			const GF61 u0 = z[k + 0], u1 = z[k + 1], u2 = z[k + 2], u3 = z[k + 3];
			z[k + 0] = u0.sqr() + u1.sqr().mul(r0); z[k + 1] = u0.mul(u1 + u1);
			z[k + 2] = u2.sqr() - u3.sqr().mul(r0); z[k + 3] = u2.mul(u3 + u3);
		}
	}

	// Radix-4, pointwise square, inverse radix-4
	void square4()
	{
		const size_t n = _vz.size();
		GF61 * const z = _vz.data();
		const GF61 * const r = _vr.data();

		for (size_t j = 0; j < n / 4; ++j)
		{
			const GF61 r0 = r[n / 4 + j];

			const size_t k = 4 * j;
			const GF61 u0 = z[k + 0], u1 = z[k + 1], u2 = z[k + 2].mul(r0), u3 = z[k + 3].mul(r0);
			const GF61 v0 = u0 + u2, v1 = u1 + u3, v2 = u0 - u2, v3 = u1 - u3;
			const GF61 s0 = v0.sqr() + v1.sqr().mul(r0), s1 = v0.mul(v1 + v1);
			const GF61 s2 = v2.sqr() - v3.sqr().mul(r0), s3 = v2.mul(v3 + v3);
			z[k + 0] = s0 + s2; z[k + 2] = (s0 - s2).mulconj(r0);
			z[k + 1] = s1 + s3; z[k + 3] = (s1 - s3).mulconj(r0);
		}
	}

	// -base < r0, r1 < base
	static GF61 reduce(int64_t & f0, int64_t & f1, const int32_t base)
	{
		int64_t l0 = f0 / base, l1 = f1 / base;
		GF61 r; r.set_int(f0 - l0 * base, f1 - l1 * base);
		f0 = l0; f1 = l1;
		return r;
	}

	// Carry propagation; output is signed -base < z[k] < base
	void carry(const bool dup)
	{
		const size_t n = _vz.size();
		GF61 * const z = _vz.data();

		const int snorm61 = _snorm61;	// 1 / (n/2)
		const int32_t base = _base;
		int64_t f0 = 0, f1 = 0;

		for (size_t k = 0; k < n; ++k)
		{
			const GF61 u = z[k] << snorm61;
			int64_t i0, i1; u.get_int(i0, i1);
			if (dup) { i0 += i0; i1 += i1; }
			f0 += i0; f1 += i1;
			z[k] = reduce(f0, f1, base);
		}

		while ((f0 != 0) || (f1 != 0))
		{
			int64_t t = f0; f0 = -f1; f1 = t;	// a_n = -a_0

			for (size_t k = 0; k < n; ++k)
			{
				int64_t i0, i1; z[k].get_int(i0, i1);
				f0 += i0; f1 += i1;
				z[k] = reduce(f0, f1, base);
				if ((f0 == 0) && (f1 == 0)) break;
			}
		}
	}

	// z = z^2 or 2z^2 modulo b^{2^n} + 1
	void step(const bool dup)
	{
		const size_t n = _vz.size();

		size_t m = n / 4, s = 1;
		for (; m > 1; m /= 4, s *= 4) forward4(m, s);
		if (m == 1) square4(); else square2x2();
		for (m = (m == 1) ? 4 : 2, s /= 4; m <= n / 4; m *= 4, s /= 4) backward4(m, s);
		carry(dup);
	}

public:
	Power(const int n)
		: _vz(size_t(1) << (n - 1)), _vr(size_t(1) << (n - 2)), _snorm61(61 - n + 2)
	{
		const size_t size = _vz.size();
		GF61 * const r = _vr.data();

		// Polynomial is z^{2^n} - (-1): init roots of -1
		for (size_t s = 1; s < size / 2; s *= 2)
		{
			const GF61 r_s = GF61::primroot_n(2 * 4 * s);
			for (size_t j = 0; j < s; ++j)
			{
				r[s + j] = r_s.pow(bitrev(j, 4 * s) + 1);
			}
		}
	}

	// 2^exponent modulo b^{2^n} + 1
	void eval(const uint32_t b, const mpz_t & exponent)
	{
		_base = b;
		GF61 * const z = _vz.data();

		z[0] = GF61(1u);
		for (size_t i = 1, size = _vz.size(); i < size; ++i) z[i] = GF61(0u);

		// Left-to-right binary exponentiation
		for (mp_bitcnt_t i = mpz_sizeinbase(exponent, 2); i > 0; --i) step(mpz_tstbit(exponent, i - 1) != 0);
	}

	// Check if output is equal to one
	bool equal_one() const
	{
		const size_t n = _vz.size();
		const GF61 * const z = _vz.data();

		std::vector<int64_t> vzi(2 * n);
		int64_t * const zi = vzi.data();

		for (size_t i = 0; i < n; ++i)
		{
			z[i].get_int(zi[i + 0 * n], zi[i + 1 * n]);
		}

		// Base-b representation: 0 <= zi[i] < b
		const int32_t base = _base;
		int64_t f;
		do
		{
			f = 0;
			for (size_t i = 0; i < 2 * n; ++i)
			{
				f += zi[i];
				int32_t r = int32_t(f % base);
				if (r < 0) r += base;
				zi[i] = r;
				f -= r;
				f /= base;
			}
			zi[0] -= f;		// a[n] = -a[0]
		} while (f != 0);

		// == 1?
		bool b = (zi[0] == 1);
		if (b) for (size_t i = 1; i < 2 * n; ++i) b &= (zi[i] == 0);
		return b;
	}
};

int main()
{
	const int n = 7;
	// GFN-7: 120, 190, 234, 506, 532, 548, 960, 1738, 1786, 2884, ...
	mpz_t exponent; mpz_init(exponent);

	Power power(n);

	// We must have 2 * 2^n * (b - 1)^2 < 2^61 - 1
	for (uint32_t b = 6, b_max = uint32_t(std::sqrt(Z61::get_p() >> (n + 1))); b <= b_max; b += 2)
	{
		if ((b != 0) && ((b & (~b + 1)) == b)) continue;	// b is a power of two => Fermat number

		// exponent is b^{2^n} + 1 - 1
		mpz_ui_pow_ui(exponent, b, static_cast<unsigned long int>(1) << n);

		power.eval(b, exponent);	// 2^exponent modulo b^{2^n} + 1

		if (power.equal_one()) std::cout << b << "^{2^" << n << "} + 1 is a probable prime" << std::endl;
	}

	mpz_clear(exponent);

	return EXIT_SUCCESS;
}
