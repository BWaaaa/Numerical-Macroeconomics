# Notebook:           Julia Function for ECON5086 Unit 4 Solving Equations
# Codes_Functions II: Solving using Newton Iteration
# Author:             2495262
# Date:               Feb 1 2021

include("Tut9_FunctionsI_Jacobianandco.jl")

## Newton Iteration: need first orther derivative
# Always, allow input x only to be arrays
function newton_iter(f::Function, x::Array{T,1},
    tol::T=1e-7, maxiter::Integer=5_000) where {T<:Real}

    count = 0
    while true
        # println(f(x))
        f′ = jacobian(f,x)
        #println("here i am", f′)
        x_new = x .- (f′ \ f(x))
        # Gausiian Elimintation must ensure f' and f are both vectors
        # println(x_new)
        if maximum(abs,x_new-x)<tol || count > maxiter
            # maximum can be broadcasted under abs operation
            break
        end
        count = count + 1
        x = copy(x_new)
    end
    print(count)

    return x,f(x), count
end

## Test: when f is good enough, then we reach the fix point quickly
# f(x)=[x[1]^2-1,x[2]^2-4]
# newton_iter(f,[1.1,1.1])
