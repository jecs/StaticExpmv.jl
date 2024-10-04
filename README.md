# StaticExpmv

[![Build Status](https://github.com/jecs/StaticExpmv.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jecs/StaticExpmv.jl/actions/workflows/CI.yml?query=branch%3Amaster)

## Description

StaticExpmv implements the operation $e^A x$ for `StaticMatrix` $A$ and `SVector` $x$ without forming $e^A$ explicitly, often resulting in dramatic speedups, especially for larger problem sizes. This implementation is based on
> Al-Mohy, A. H., & Higham, N. J. (2011). Computing the action of the matrix exponential, with an application to exponential integrators. SIAM journal on scientific computing, 33(2), 488-511.

## Usage
```
using StaticExpmv,StaticArrays
import LinearAlgebra
A = I+SMatrix{10,10}(randn(10,10))
x = SVector{10}(randn(10))
y = expmv(A,x)
```
