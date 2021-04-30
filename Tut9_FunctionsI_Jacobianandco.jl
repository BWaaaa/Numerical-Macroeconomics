# Notebook:          III Macro_6 Computational Macro Unit 3
# Codes_Functions I: Jacobian and derivatives
# Author:            2495262
# Date:              Jan 6 2021

function gradient(f::Function,x::Array{T,1}) where {T<:Real}

    m = length(x)

	h  = (eps(Float64)/2)^(1/3)*maximum(abs,[x;1.0])
    dh = Matrix(1.0I,m,m)*h

	deriv = Array{Float64}(undef,1,m)

    for i = 1:m
	    f1 = f(x.+dh[:,i])
	    f2 = f(x.-dh[:,i])
        deriv[1,i] = (f1-f2)/(2.0*h)
    end

    return deriv

end


function hessian(f::Function,x::Array{T,1}) where {T<:Real}

    m = length(x)

	h = (4*eps(Float64))^(1/3)*maximum(abs,[x;1.0])
	dh = Matrix(1.0I,m,m)*h

    hess = Array{Float64}(undef,m,m)

    for i = 1:m
        for j = i:m
	        f1 = f(x.+dh[:,i].+dh[:,j])
	        f2 = f(x.-dh[:,i].+dh[:,j])
	        f3 = f(x.+dh[:,i].-dh[:,j])
	        f4 = f(x.-dh[:,i].-dh[:,j])
	        hess[i:i,j:j] .= (f1-f2-f3+f4)/(4*h^2)
            if i != j
                hess[j,i] = hess[i,j]
            end
        end
    end

    return hess

end


function jacobian(f::Function,x::Array{T,1}) where {T <:Real}

    n = length(f(x))
    m = length(x)

	h  = (eps(Float64)/2)^(1/3)*maximum(abs,[x;1.0])
    dh = Matrix(1.0I,m,m)*h

	deriv = Array{Float64}(undef,n,m)

    for i = 1:m
	    f1 = f(x.+dh[:,i])
	    f2 = f(x.-dh[:,i])
        deriv[:,i] .= (f1-f2)/(2.0*h)
    end

    return deriv

end
