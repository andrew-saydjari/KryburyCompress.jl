"""
    kryburyCompress(M::LowRankMultMat,dM::LowRankDiagMat,nfeat::Int;nsub=2,tol=1e-12) -> W

Given a matrix M, finds an approximation M ≈ A + V*V' using Krylov methods to approximate M by `nsub` eigenvectors and assign the residual diagonal to A. To aid in speed, M is passed as a two objects of type `LowRankMultMat` and `LowRankDiagMat.`

# Arguments:
- `M`: LowRankMultMat representation of M
- `dM`: LowRankDiagMat representation of M
- `nfeat`: linear dimension of M (one side length)

# Keywords:
- `nsub`: number of eigenvectors in approximation of M
- `tol`: tolerance of Krylov approximation

# Outputs:
- `W`: DiagWoodbury object that approximates M
"""
function kryburyCompress(M::LowRankMultMat,dM::LowRankDiagMat,nfeat::Int;nsub=2,tol=1e-12)
    λ, V, info = eigsolve(M,
        nfeat,
        nsub,
        krylovdim=2*nsub,
        ishermitian=true,
        issymmetric=true,
        :LM,
        tol=tol
    );
    Vmat = hcat(V[1:nsub]...)*Diagonal(sqrt.(λ[1:nsub]))
    true_diag = diag(dM)
    A = abs.(true_diag.-dropdims(sum(Vmat.^2,dims=2),dims=2))
    return DiagWoodbury(A,Vmat,true_diag)
end

"""
    save(W::DiagWoodbury,fname::String)

Writes DiagWoodbury type to a FITS file.

# Arguments:
- `W`: DiagWoodbury object to be saved
- `fname`: name of FITS file to be written
"""
function save(W::DiagWoodbury,fname::String)
    if last(fname,5) == ".fits"
        fname = chop(fname,tail=5)
    end
    f = FITS(fname*".fits", "w");
    
    hdrA = FITSIO.default_header(W.A)
    hdrA["NAME"] = "A"
    set_comment!(hdrA, "NAME", "diagonal component of Woodbury form")
    write(f,W.A,header=hdrA,name="A")
    
    hdrV = FITSIO.default_header(W.V)
    hdrV["NAME"] = "V"
    set_comment!(hdrV, "NAME", "low-rank component of Woodbury form")
    write(f,W.V,header=hdrV,name="V")
    
    hdrD = FITSIO.default_header(W.D)
    hdrD["NAME"] = "D"
    set_comment!(hdrD, "NAME", "best diagonal approximation")
    write(f,W.D,header=hdrD,name="D")
    close(f)
end

"""
    read_krybury(fname::String) -> W

Reads a DiagWoodbury type from a FITS file.

# Arguments:
- `fname`: name of FITS file to be read

# Outputs:
- `W`: DiagWoodbury object from file
"""
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
