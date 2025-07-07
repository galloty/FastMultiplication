# Fast Multiplication

## About

Fast integer multiplication is an algorithm able to perform the multiplication of two *n*-digit integers using *O*(*n*&nbsp;log&nbsp;*⁡n*) operations.

Modular multiplication of polynomials is the basis of this method.  
An intricate method is often used: the product can be rewritten as a cyclic convolution and is evaluated using Fast Fourier Transforms. It is unnecessary and the purpose of this project is to demystify the fast multiplication for large integers.

Fast modular multiplication of polynomials is a divide-and-conquer algorithm, where factorization and the Chinese remainder theorem are applied recursively. This method is described here: starting from a simple implementation of the mathematical relation, a list of elementary programs illustrates how it can be extended and improved.

## [*fastMul_rec.cpp*](fastMul_rec.cpp)

The recursion is applied to compute *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>8</sup>&nbsp;&minus;&nbsp;7, where *P* is a polynomial of degree 8&minus;1. The field is the complex numbers.  
*P*(*x*) mod&nbsp;*x*<sup>*m*</sup>&nbsp;&minus;&nbsp;*r*<sup>*m*</sup> and *P*(*x*) mod&nbsp;*x*<sup>*m*</sup>&nbsp;+&nbsp;*r*<sup>*m*</sup> are calculated from *P*(*x*) mod&nbsp;*x*<sup>2*m*</sup>&nbsp;&minus;&nbsp;*r*<sup>2*m*</sup>. The last step is Lagrange interpolation and we have *P*(*x*) mod&nbsp;*x*&nbsp;&minus;&nbsp;*r*&nbsp;= *P*(*r*). *P*(*r*)<sup>2</sup> is a scalar. Finally, *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>2*m*</sup>&nbsp;&minus;&nbsp;*r*<sup>2*m*</sup> is calculated from *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>*m*</sup>&nbsp;&minus;&nbsp;*r*<sup>*m*</sup> and *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>*m*</sup>&nbsp;+&nbsp;*r*<sup>*m*</sup> with the Chinese remainder theorem for polynomials.

## [*fastMul_rec_235.cpp*](fastMul_rec_235.cpp)

The algorithm is extended to *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>2<sup>2</sup>&nbsp;&middot;&nbsp;3<sup>2</sup>&nbsp;&middot;&nbsp;5<sup>2</sup></sup>&nbsp;&minus;&nbsp;7. We have *x*<sup>3*m*</sup>&nbsp;&minus;&nbsp;*r*<sup>3*m*</sup>&nbsp;= (x<sup>*m*</sup>&nbsp;&minus;&nbsp;*r*<sup>*m*</sup>)&nbsp;&middot;&nbsp;(x<sup>*m*</sup>&nbsp;&minus;&nbsp;*jr<sup>*m*</sup>*)&nbsp;&middot;&nbsp;(x<sup>*m*</sup>&nbsp;&minus;&nbsp;*j*<sup>2</sup>*r*<sup>*m*</sup>)&nbsp;&middot;&nbsp;, where *j* is a primitive root of *x*<sup>3</sup>&nbsp;&minus;&nbsp;1</sup>. The equivalent relation is applied to *x*<sup>5*m*</sup>&nbsp;&minus;&nbsp;*r*<sup>5*m*</sup>.

## [*fastMul_rec_GF.cpp*](fastMul_rec_GF.cpp)

Roots are calculated over a prime finite field of order *p*. There are two conditions:  
 - If *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>*n*</sup>&nbsp;&minus;&nbsp;*r* is evaluated then a *n*<sup>th</sup> root of *r* must be an element in Z/*p*Z.  
 - If the coefficients of *P*(*x*) are non-negative integers smaller than *b* (*P*(*b*) is the base-*b* representation of an integer) then *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>*n*</sup>&nbsp;&minus;&nbsp;*r* over the natural numbers can be retrieved from the result over Z/*p*Z if *r*&nbsp;> 0 and *p*&nbsp;> *n*&nbsp;&middot;&nbsp;(*b*&nbsp;&minus;&nbsp;1)<sup>2</sup>. If *r*&nbsp;< 0 then outputs can be negative integers: the condition is *p*&nbsp;> 2&nbsp;*n*&nbsp;&middot;&nbsp;(*b*&nbsp;&minus;&nbsp;1)<sup>2</sup> to determine the sign of each coefficient.  

Here *n*&nbsp;= 900, *r*&nbsp;= 7 and *b*&nbsp;= 1000. *p*&nbsp;= 913262401 is chosen: we have 243771734<sup>*n*</sup>&nbsp;= 7 (mod&nbsp;*p*) and *p*&nbsp;> *n*&nbsp;&middot;&nbsp;(*b*&nbsp;&minus;&nbsp;1)<sup>2</sup>&nbsp;= 898200900. Rather than evaluating the square, cube or 5<sup>th</sup> root of an element, the *n*<sup>th</sup> root of unity and the *n*<sup>th</sup> root of 7 are computed and roots needed for the algorithm are of the form 1<sup>*u*/*n*</sup>&nbsp;&middot;&nbsp;7<sup>*v*/*n*</sup>.

## [*fastMul_rec_twisted.cpp*](fastMul_rec_twisted.cpp)

The recursion is based on the relation *x*<sup>*am*</sup>&nbsp;&minus;&nbsp;*r*<sup>*am*</sup>&nbsp;= Prod<sub>0&le;*k*<*a*</sub>&nbsp;(*x*<sup>*m*</sup>&nbsp;&minus;&nbsp;&alpha;<sup>*k*</sup>*r*<sup>*m*</sup>), where *a* is prime and &alpha; is a primitive root of *x*<sup>*a*</sup>&nbsp;&minus;&nbsp;1</sup>.  
Let *X*&nbsp;= *x*<sup>*m*</sup> and *R*&nbsp;= *r*<sup>*m*</sup>. If *X*&nbsp;&minus;&nbsp;&alpha;<sup>*k*</sup>*R* is *twisted* into *Y*<sub>*k*</sub>&nbsp;&minus;&nbsp;*R* then at each step all polynomial moduli are identical.  
*P*(*X*)&nbsp;mod&nbsp;*X*&nbsp;&minus;&nbsp;&alpha;<sup>*k*</sup>*R*&nbsp;= *P*(*X*) mod&nbsp;&alpha;<sup>&minus;*k*</sup>*X*&nbsp;&minus;&nbsp;*R*. Let *Y*&nbsp;= &alpha;<sup>&minus;*k*</sup>*X*, we have *y*&nbsp;= &alpha;<sup>&minus;*k*/*m*</sup>*x*.  
If *P*(*x*)&nbsp;= Sum<sub>0&le;*i*<*m*</sub>&nbsp;c<sub>*i*</sub>&nbsp;*x*<sup>*i*</sup> then *P*(*y*)&nbsp;= Sum<sub>0&le;*i*<*m*</sub>&nbsp;c<sub>*i*</sub>&nbsp;&alpha;<sup>*ik*/*m*</sup>&nbsp;*y*<sup>*i*</sup>. The coefficients of *P*(*y*) are *c*'<sub>*i*</sub>&nbsp;= &alpha;<sup>*ik*/*m*</sup>&nbsp;*c*<sub>*i*</sub>. After the recursion, the coefficients of *P*(*y*) are *untwisted* into *P*(*x*) with *c*<sub>*i*</sub>&nbsp;= &alpha;<sup>&minus;*ik*/*m*</sup>&nbsp;*c*'<sub>*i*</sub>.  

Twisting is the generalization of weighted transforms, where weighting is applied at each step.  

With the previous version of the fast multiplication, *n*/2 multiplications by 1<sup>*u*/*n*</sup>&nbsp;&middot;&nbsp;*r*<sup>*v*/*n*</sup> are computed at each step. With the twisted form, *n*/2 multiplications by *r*<sup>*v*/*n*</sup> and *n*/2 multiplications by 1<sup>*u*'/*n*</sup> are needed. The twisted transform is as efficient as the non-twisted form only if *r*&nbsp;= 1.

Twisted FFT was defined in Daniel J. Bernstein, [Multidigit multiplication for mathematicians](https://cr.yp.to/papers/m3-20010811-retypeset-20220327.pdf), 2001. It was restricted to *r*&nbsp;= 1 and *n* is a power of two.

## [*fastMul_rec_FFT.cpp*](fastMul_rec_FFT.cpp)

The twisted fast multiplication is written for *r*&nbsp;= 1. Roots of *r* are equal to one and it can be seen that the twist factors &alpha;<sup>*ik*/*m*</sup> are the twiddle factors of a Fast Fourier Transform.  

The *transform* that splits *P*(*x*) mod&nbsp;*x*<sup>*n*</sup>&nbsp;&minus;&nbsp;1 into *P*(*x*) mod&nbsp;*x*&nbsp;&minus;&nbsp;&omega;<sup>*i*</sup>, where &omega; is a primitive *n*<sup>th</sup> root of unity and 0&le;*i*<*n* is Gentleman-Sande FFT recursion. The inverse transform that merges *P*(*x*)<sup>2</sup> mod&nbsp;*x*&nbsp;&minus;&nbsp;&omega;<sup>*i*</sup> into *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>*n*</sup>&nbsp;&minus;&nbsp;1 is Cooley-Tukey FFT recursion. Because Gentleman-Sande and Cooley-Tukey diagrams are symmetric, data ordering of the *output* is irrelevant, even if an iterative FFT algorithm is implemented.

The non-twisted fast multiplication can also be used for the computation of *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>*n*</sup>&nbsp;&minus;&nbsp;1. It is a different algorithm.  
The forward transform was found by Georg Bruun, ["z-Transform DFT filters and FFTs](https://backend.orbit.dtu.dk/ws/files/4658740/Bruun.pdf), 1978, IEEE Transactions on Acoustics, Speech, and Signal Processing, 26 (1): 56-63. Here, the Discrete Fourier Transform is written as a Z-transform filter. A Discrete Fourier Transform is evaluated with *n* roots of unity along the z-domain's unit circle. But this restriction does not apply to Z-transforms. If the FFT based multiplication is limited to the polynomial *x*<sup>n</sup>&nbsp;&minus;&nbsp;1 (cyclic convolution), Bruun's method can be applied to any polynomial.  

The two algorithms are implemented: because *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>*n*</sup>&nbsp;&minus;&nbsp;1 is calculated, both are almost equivalent. But the method based on polynomial factorization can compute *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>*n*</sup>&nbsp;&minus;&nbsp;*r* or *P*(*x*)<sup>2</sup> mod&nbsp;&Phi;<sub>*m*</sub>(*x*<sup>*n*</sup>), where &Phi;<sub>*m*</sub> is the *m*<sup>th</sup> cyclotomic polynomial, without an extra "weighted" step.

## [*fastMul_GF.cpp*](fastMul_GF.cpp)

Fast multiplication based on a recursive polynomial factorization is an in-place algorithm. The roots needed at each step can be precomputed. An iterative (non-recursive) version of [fastMul_rec_GF.cpp](fastMul_rec_GF.cpp) is implemented here. Note that roots are stored in 2-3-5-reversed order: it is a generalization of bit-reversal permutation based on decomposing *n* into its prime factors. Rather than a division by 2, 3 or 5 at each step of the reverse process, a single division by *n* is computed at the end of the algorithm.  

The number of operations is equal to a FFT based multiplication, where the polynomial of the convolution is *x*<sup>900</sup>&nbsp;&minus;&nbsp;1. But here the polynomial is *x*<sup>900</sup>&nbsp;&minus;&nbsp;7.
