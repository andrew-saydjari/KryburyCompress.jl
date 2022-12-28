# KryburyCompress

[![][action-img]][action-url]
[![][codecov-img]][codecov-url]

The purpose of this package is to enable the compression of possibly high dimensional covariance matrices that are well represented by low-rank approximations. In general, when covariances in data are large, one must store either an N x N covariance matrix or, one makes the gross approximation that the covariance matrix is diagonal. This package allows interpolation between those two limits, finding a fast, low rank approximation to the covariance matrix which is exported as a type with fast matrix-vector multiplication. We provide read/write methods to store this object in a FITS file, to allow multilingual compatability.

## Installation

Currently, installation is directly from the GitHub

```julia
import Pkg
Pkg.add(url="https://github.com/andrew-saydjari/KryburyCompress.jl")
```

## Compression

Given a matrix $M$, we find an approximation $M \approxeq A + V * V'$ where $A$ is `Diagonal` and $V$ is a skinny matrix, containing the largest eigenvectors of $M$. The largest eigenvectors are found in-practice using Krylov methods from [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl). Then, $A$ is simply the residual diagonal between $M$ and $V * V'$. This compression can be massively accelerated when $M$ is itself composed of products of low-rank matrices represented via [LowRankOps.jl](https://github.com/andrew-saydjari/LowRankOps.jl).

### Example Usage:

We use a `LowRankMultMat` and `LowRankDiagMat` here to speed up the compression. In this example, $M = (U * U') - (U * U') * D * (U * U')$ and we want to find a low-rank approximation to $M$.

```julia
using KryburyCompress
using LinearAlgebra, KrylovKit, FITSIO, LowRankOps
using Random: seed!

seed!(123)
n_big = 100
n_lit = 2

D = 1e-4*abs.(randn(n_big))
U = randn(n_big,n_lit)

function Cii_precomp_mult(matList)
    return []
end

function Cii_fxn_mult(matList,precompList,x)
    Ctotinv = matList[1]
    Vi = matList[2]
    arg2 = Vi*(Vi'*x)
    return arg2 - Vi*(Vi'*(Ctotinv*arg2))
end

function Cii_precomp_diag(matList)
    Ctotinv = matList[1]
    Vi = matList[2]
    return [Vi'*(Ctotinv*Vi)]
end

function Cii_diag_map(matList,precompList)
    Vi = matList[2]
    arg1 = precompList[1]
    return dropdims(sum(Vi.^2,dims=2),dims=2).-dropdims(sum(Vi'.*(arg1*Vi'),dims=1),dims=1)
end

BMat = LowRankMultMat([Diagonal(D),U],Cii_precomp_mult,Cii_fxn_mult);
BMatDiag = LowRankDiagMat([Diagonal(D),U],Cii_precomp_diag,Cii_diag_map);

W = kryburyCompress(BMat,BMatDiag,n_big)
```

## I/O

The `DiagWoodbury` object resulting from the compression can be written and read to a FITS file using 

```julia
fname = "example.fits"
save(W,fname)
W = read_krybury(fname)
```

Each of the fields from the `DiagWoodbury` object are saved as different extensions in the FITS file.

<!-- URLS -->
[action-img]: https://github.com/andrew-saydjari/KryburyCompress.jl/workflows/CI/badge.svg
[action-url]: https://github.com/andrew-saydjari/KryburyCompress.jl/actions

[codecov-img]: https://codecov.io/github/andrew-saydjari/KryburyCompress.jl/coverage.svg?branch=main
[codecov-url]: https://codecov.io/github/andrew-saydjari/KryburyCompress.jl?branch=main
