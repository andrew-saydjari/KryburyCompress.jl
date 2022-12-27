@testset "diagWoodbury.jl" begin
    seed!(123)
    n_big = 100
    n_lit = 2
    
    A = randn(n_big)
    V1 = randn(n_big,n_lit)
    x = randn(n_big)
    X = randn(n_big,n_big)
    
    W = DiagWoodbury(A,V1)
    W1 = DiagWoodbury(A,V1,diag(Diagonal(A)+V1*V1'))
    FM = Diagonal(A) + V1*V1'
    
    @test W ≈ W1
    @test W + W1 == DiagWoodbury(2*A,[V1 V1])
    @test (2*W)*x ≈ 2*(W*x)
    @test W*X ≈ FM*X
    @test W*x ≈ FM*x
end