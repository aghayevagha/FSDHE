# Fourier Shape Descriptors with Homomorphic Encryption

The repository contains the implementation of an efficient shape matching algorithm in homomorphic domain, with OpenFHE library.



## Overview
The algorithm takes the boundary points (2D) of two shapes, and compares their closeness with frequency analysis. It converts two-dimensional points to 
Fourier Shape descriptors, and compute the Euclidian difference between them, if the difference is smaller than one, we can assume they are close, otherwise they
have different shapes. We have provided two examples, and detailed explanation can be found on our paper at https://eprint.iacr.org/2025/1942




## Usage
You should install OpenFHE library for C++ language, and the easiest way of running this code is to clone the repo to `src\pke` and run following:
```
make 
cmake ..
./bin/examples/FSDHE/FSDHE.cpp
```
You can refer to https://jdumezy.com/blog/openfhe-crash-course/ for detailed explanation on how to build the library.
