# Fast Multiplication

## About

Fast integer multiplication is an algorithm able to perform the multiplication of two *n*-digit integers using *O*(*n*&nbsp;log&nbsp;*⁡n*) operations.

Modular multiplication of polynomials is the basis of this method.  
A complex approach is often used: the product can be rewritten as a cyclic convolution and is evaluated using Fast Fourier Transforms. It is unnecessary and the purpose of this project is to demystify the fast multiplication for large integers.

Fast modular multiplication of polynomials is a divide-and-conquer algorithm, where the Chinese remainder theorem is applied recursively. This method is described here: starting from a simple implementation of the mathematical relation, a list of elementary programs illustrates how it can be improved.

## [*fastMul_rec.cpp*](fastMul_rec.cpp)

The recursion is applied to compute *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>8</sup>&nbsp;&minus;&nbsp;7, where *P* is a polynomial of degree 8&minus;1. The field is the complex numbers. The Chinese remainder theorem for polynomials is evaluated: at each step, *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>2*m*</sup>&nbsp;&minus;&nbsp;*r*<sup>2*m*</sup> is calculated from *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>*m*</sup>&nbsp;&minus;&nbsp;*r*<sup>*m*</sup> and *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>*m*</sup>&nbsp;+&nbsp;*r*<sup>*m*</sup>. The last step is Lagrange interpolation and we have *P*(*x*) mod&nbsp;*x*&nbsp;&minus;&nbsp;*r* = *P*(*r*).

## [*fastMul_rec235.cpp*](fastMul_rec235.cpp)

The algorithm is extended to *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>2<sup>2</sup>&nbsp;&middot;&nbsp;3<sup>2</sup>&nbsp;&middot;&nbsp;5<sup>2</sup></sup>&nbsp;&minus;&nbsp;7. We have *x*<sup>3*m*</sup>&nbsp;&minus;&nbsp;*r*<sup>3*m*</sup> = (x&nbsp;&minus;&nbsp;*r*)&nbsp;&middot;&nbsp;(x&nbsp;&minus;&nbsp;*jr*)&nbsp;&middot;&nbsp;(x&nbsp;&minus;&nbsp;*j*<sup>2</sup>*r*)&nbsp;&middot;&nbsp;, where *j* is a primitive root of *x*<sup>3</sup>&nbsp;&minus;&nbsp;1</sup>. The equivalent relation is applied to *x*<sup>5*m*</sup>&nbsp;&minus;&nbsp;*r*<sup>5*m*</sup>.

## [*fastMul_rec_GF.cpp*](fastMul_rec_GF.cpp)

Roots are calculated over a prime finite field of order *p*. There are two conditions:  
 - If *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>*n*</sup>&nbsp;&minus;&nbsp;*r* is evaluated then a *n*<sup>th</sup> root of *r* must be an element in Z/*p*Z.  
 - If the coefficients of *P*(*x*) are non-negative integers smaller than *b* (*P*(*x*) is the base-*b* representation of an integer) then *P*(*x*)<sup>2</sup> mod&nbsp;*x*<sup>*n*</sup>&nbsp;&minus;&nbsp;*r* over the natural numbers can be retrieved from the result over Z/*p*Z if *r*&nbsp;> 0 and *p*&nbsp;> *n*&nbsp;&middot;&nbsp;(*b*&nbsp;&minus;&nbsp;1)<sup>2</sup>. If *r*&nbsp;< 0 then outputs can be negative integers: the condition is *p*&nbsp;> 2&nbsp;*n*&nbsp;&middot;&nbsp;(*b*&nbsp;&minus;&nbsp;1)<sup>2</sup> to determine the sign of each coefficient.  

Here *n*&nbsp;= 900, *r*&nbsp;= 7 and *b*&nbsp;= 1000. *p*&nbsp;= 913262401 is chosen: we have 243771734<sup>*n*</sup>&nbsp;= 7 (mod&nbsp;*p*) and *p*&nbsp;> *n*&nbsp;&middot;&nbsp;(*b*&nbsp;&minus;&nbsp;1)<sup>2</sup>&nbsp;= 898200900. Rather than evaluating the square, cube or 5<sup>th</sup> root of an element, the *n*<sup>th</sup> root of unity and the *n*<sup>th</sup> root of 7 are computed and roots needed for the algorithm are of the form 1<sup>*u*/*n*</sup>&nbsp;&middot;&nbsp;7<sup>*v*/*n*</sup>.

