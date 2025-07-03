# Fast Multiplication

## About

Fast integer multiplication is an algorithm able to perform the multiplication of two *n*-digit integers using *O*(*n* log *⁡n*) operations.

Modular multiplication of polynomials is the basis of this method.  
A complex approach is often used: the product can be rewritten as a cyclic convolution and is evaluated using Fast Fourier Transforms. It is unnecessary and the purpose of this project is to demystify the fast multiplication for large integers.

Fast modular multiplication of polynomials is a divide-and-conquer algorithm, where the Chinese remainder theorem is applied recursively. This method is described here: starting from a simple implementation of the mathematical relation, a list of elementary programs illustrates how it can be improved.

## *fastMul_rec.cpp*

The recursion is applied to compute *P*(*x*)<sup>2<sup> mod *x*<sup>8<sup> - 7, where *P* is a polynomial of degree 8-1. The field is the complex numbers. The Chinese remainder theorem for polynomials is evaluated: at each step, *P*(*x*)<sup>2<sup> mod *x*<sup>2*m*<sup> - *r*<sup>2*m*<sup> is calculated from *P*(*x*)<sup>2<sup> mod *x*<sup>*m*<sup> - *r*<sup>*m*<sup> and *P*(*x*)<sup>2<sup> mod *x*<sup>*m*<sup> + *r*<sup>*m*<sup>. The last step is Lagrange interpolation and we have *P*(*x*) mod *x* - *r*<sup> = *P*(*r*).

## *fastMul_rec235.cpp*

