struct DiagWoodbury{T<:Vector,U<:Matrix,S<:Vector}
    A::T
    V::U
    D::S
end

DiagWoodbury(A::Vector,V::Matrix) = DiagWoodbury(A,V,A+dropdims(sum(V.^2,dims=2),dims=2));

==(x::DiagWoodbury,y::DiagWoodbury) = ((x.A==y.A) & (x.V==y.V)) & (x.D==y.D)

+(W::DiagWoodbury, X::DiagWoodbury) = DiagWoodbury(W.A + X.A, [W.V X.V])

*(α::Real, W::DiagWoodbury)= DiagWoodbury(α*W.A, sqrt(α)*W.V)

*(W::DiagWoodbury, B::AbstractMatrix)=Diagonal(W.A)*B + W.V*(W.V'*B)

*(W::DiagWoodbury, B::Vector)=Diagonal(W.A)*B + W.V*(W.V'*B)