module KryburyCompress

    using LinearAlgebra, KrylovKit, FITSIO, LowRankOps
    import Base: *, +, ==, ≈

    export DiagWoodbury, save, read_krybury, kryburyCompress

    include("diagWoodbury.jl")
    include("compress.jl")

end
