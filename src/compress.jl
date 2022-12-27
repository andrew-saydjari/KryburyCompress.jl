function kryburyCompress(M::LowRankMultMat,dM::LowRankDiagMat,nfeat;nsub=2,tol=1e-12)
    λ, V, info = eigsolve(M,
        nfeat,
        nsub,
        krylovdim=2*nsub,
        ishermitian=true,
        issymmetric=true,
        :LM);
    Vmat = hcat(V[1:nsub]...)*Diagonal(sqrt.(λ[1:nsub]))
    true_diag = diag(dM)
    A = abs.(true_diag.-dropdims(sum(Vmat.^2,dims=2),dims=2))
    return DiagWoodbury(A,Vmat,true_diag)
end

function save(W::DiagWoodbury,fname::String)
    if last(fname,5) == ".fits"
        fname = chop(fname,tail=5)
    end
    f = FITS(fname*".fits", "w");
    
    hdrA = FITSIO.default_header(W.A)
    hdrA["NAME"] = "A"
    set_comment!(hdrA, "NAME", "diagonal component of Woodbury form")
    write(f,W.A,header=hdrA,name="A")
    
    hdrV = FITSIO.default_header(out.V)
    hdrV["NAME"] = "V"
    set_comment!(hdrV, "NAME", "low-rank component of Woodbury form")
    write(f,W.V,header=hdrV,name="V")
    
    hdrD = FITSIO.default_header(out.D)
    hdrD["NAME"] = "D"
    set_comment!(hdrD, "NAME", "best diagonal approximation")
    write(f,W.D,header=hdrD,name="D")
    close(f)
end

function read_krybury(fname::String)
    if last(fname,5) == ".fits"
        fname = chop(fname,tail=5)
    end
    f = FITS(fname*".fits", "r+");
    A = read(f["A"])
    V = read(f["V"])
    D = read(f["D"])
    close(f)
    return DiagWoodbury(A,V,D)
end