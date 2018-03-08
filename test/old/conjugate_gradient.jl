if VERSION >= v"0.7-"
    using Random
end

@testset "conjug_grad" begin

srand(42)   # always use same random numbers
n = 40
A = rand(n, n)
A = 0.5*(A + A')
b = rand(n)
y = rand(n)
ɛ = 1e-5

x = conjugate_gradient(A, b, y; ɛ=ɛ, max_iter=1000)

@test norm2(A*x - b) <= ɛ * norm(A, Inf) * norm(b, Inf)

end
 
