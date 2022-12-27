module KryburyCompress

    using LinearAlgebra, FITSIO
    import Base: *, +, ==

    export DiagWoodbury, save, read_krybury, kryburyCompress

    include("diagWoodbury.jl")
    include("compress.jl")

end
