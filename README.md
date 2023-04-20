# Fast Multiplication

## About

Fast integer multiplication is an algorithm able to perform the multiplication of two *n*-digit integers using *O*(*n* log *⁡n*) operations.

Modular multiplication of polynomials is the basis of this method.  
A complex approach is often used: the product can be rewritten as a cyclic convolution and is evaluated using a Fast Fourier Transform. It is unnecessary and the purpose of this project is to demystify fast multiplication algorithm for large integers.

Fast modular multiplication of polynomials is a divide-and-conquer algorithm, where the Chinese remainder theorem is applied recursively. This method is described here: starting from a simple implementation of the mathematical relation, a list of elementary programs illustrate how it can be improved.
