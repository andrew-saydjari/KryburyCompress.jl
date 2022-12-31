@testset "diagWoodbury.jl" begin
    seed!(123)
    n_big = 100
    n_lit = 2
    
    A = 1e-4*abs.(randn(n_big))
    V = randn(n_big,n_lit)
    x = randn(n_big)
    X = randn(n_big,n_big)
    
    W = DiagWoodbury(A,V)
    W1 = DiagWoodbury(A,V,diag(Diagonal(A)+V*V'))
    FM = Diagonal(A) + V*V'
    
    @test W ≈ W1
    @test W + W1 == DiagWoodbury(2*A,[V V])
    @test (2*W)*x ≈ 2*(W*x)
    @test W*X ≈ FM*X
    @test W*x ≈ FM*x
    @test inv(FM)*x ≈ W\x
    @test inv(FM)*X ≈ W\X
    @test Matrix(W) == FM
    @test W' == W
    @test diag(W) == W.D
end