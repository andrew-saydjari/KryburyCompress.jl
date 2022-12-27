@testset "compress.jl" begin
    seed!(123)
    n_big = 100
    n_lit = 2
    
    A = 1e-4*abs.(randn(n_big))
    V = randn(n_big,n_lit)
    W = DiagWoodbury(A,V)
    
    fname = tempname()
    save(W,fname)
    W1 = read_krybury(fname)
    @test W == W1
    
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
    
    BMat = LowRankMultMat([Diagonal(A),V],Cii_precomp_mult,Cii_fxn_mult);
    BMatDiag = LowRankDiagMat([Diagonal(A),V],Cii_precomp_diag,Cii_diag_map);
    
    out = kryburyCompress(BMat,BMatDiag,n_big)
    
    fname = tempname()
    save(out,fname)
    W2 = read_krybury(fname)
    @test out == W2
    
end