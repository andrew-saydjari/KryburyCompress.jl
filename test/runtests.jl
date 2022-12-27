using KryburyCompress
using LinearAlgebra, KrylovKit, FITSIO, LowRankOps, Test
using Random: seed!

include("diagWoodbury.jl")
include("compress.jl")