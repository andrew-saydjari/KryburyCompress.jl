@testset "compress.jl" begin
    seed!(123)
    n_big = 100
    n_lit = 2
    
    A = 1e-4*abs.(randn(n_big))
    V = randn(n_big,n_lit)
    W = DiagWoodbury(A,V)
    
    fname = tempname()*".fits"
    save(W,fname)
    W1 = read_krybury(fname)
    @test W == W1
    
    function dummy_precomp_mult(matList)
        return []
    end

    function dummy_fxn_mult(matList,precompList,x)
        A = matList[1]
        Vi = matList[2]
        return A*x + Vi*(Vi'*x)
    end

    function dummy_precomp_diag(matList)
        return []
    end

    function dummy_diag_map(matList,precompList)
        A = matList[1]
        Vi = matList[2]
        return diag(A) + dropdims(sum(Vi.^2,dims=2),dims=2)
    end
    
    BMat = LowRankMultMat([Diagonal(A),V],dummy_precomp_mult,dummy_fxn_mult);
    BMatDiag = LowRankDiagMat([Diagonal(A),V],dummy_precomp_diag,dummy_diag_map);
    
    out = kryburyCompress(BMat,BMatDiag,n_big)
    
    fname = tempname()
    save(out,fname)
    W2 = read_krybury(fname)
    @test out == W2
    
end
