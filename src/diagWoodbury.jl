"""
This is a `DiagWoodbury` object.
It can be constructed from A and V alone.
Implemented methods include `==`, `≈`, `+`, `*`, `\`, `Matrix()`, `adjoint()`, `diag()`.

# Fields
- A: diagonal component of the decomposition of M
- V: second letter of the English alphabet
- D: best diagonal approximation the matrix M
- invPrecomp: precomputation used for `\` implementation
"""
struct DiagWoodbury{T<:Vector,U<:Matrix,S<:Vector,X<:Matrix}
    A::T
    V::U
    D::S
    invPrecomp::X
end

function DiagWoodbury(A::AbstractVector,V::Matrix,D::AbstractVector)
    Ainv = Diagonal(1 ./A)
    AinvV = Ainv*V
    return DiagWoodbury(A,
        V,
        D,
        AinvV*inv(I + V'*AinvV)
        );
end

function DiagWoodbury(A::AbstractVector,V::Matrix)
    Ainv = Diagonal(1 ./A)
    AinvV = Ainv*V
    return DiagWoodbury(A,
        V,
        A+dropdims(sum(V.^2,dims=2),dims=2),
        AinvV*inv(I + V'*AinvV)
        );
end

==(x::DiagWoodbury,y::DiagWoodbury) = ((x.A==y.A) & (x.V==y.V)) & (x.D==y.D)
≈(x::DiagWoodbury,y::DiagWoodbury) = ((x.A≈y.A) & (x.V≈y.V)) & (x.D≈y.D)

+(W::DiagWoodbury, X::DiagWoodbury) = DiagWoodbury(W.A + X.A, [W.V X.V])

*(α::Real, W::DiagWoodbury)= DiagWoodbury(α*W.A, sqrt(α)*W.V)

*(W::DiagWoodbury, B::AbstractMatrix)=Diagonal(W.A)*B + W.V*(W.V'*B)
*(W::DiagWoodbury, B::AbstractVector)=Diagonal(W.A)*B + W.V*(W.V'*B)

function \(W::DiagWoodbury, X::AbstractVector)
    Ainv = Diagonal(1 ./W.A)
    return Ainv*(X - W.V*(W.invPrecomp'*X))
end

function \(W::DiagWoodbury, X::AbstractMatrix)
    Ainv = Diagonal(1 ./W.A)
    return Ainv*(X - W.V*(W.invPrecomp'*X))
end

Matrix(W::DiagWoodbury)= Diagonal(W.A) + W.V*W.V'

adjoint(W::DiagWoodbury)= W

diag(W::DiagWoodbury)= W.D