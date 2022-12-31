module KryburyCompress

    using LinearAlgebra, KrylovKit, FITSIO, LowRankOps
    import Base: *, +, ==, ≈, \
    import LinearAlgebra: Matrix, adjoint, diag

    export DiagWoodbury, save, read_krybury, kryburyCompress

    include("diagWoodbury.jl")
    include("compress.jl")

end
