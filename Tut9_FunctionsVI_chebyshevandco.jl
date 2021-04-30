# Notebook:           Julia Function for ECON5086 Unit 4 Solving Equations
# Codes_Functions IV: Lecture 7 Integration Question 5
# Author:             2495262
# Date:               Mar 16 2021

include("Tut9_FunctionsIII_chebyshevAppx.jl")


function chebyshev_weights(y::Array{Float64,2},nodes1::Array{Float64,1},nodes2::Array{Float64,1},order::Array{Int64,1},domain::Array{Float64,2})

  nodes1 = 2.0*(nodes1.-domain[2,1])/(domain[1,1]-domain[2,1]).-1.0
  nodes2 = 2.0*(nodes2.-domain[2,2])/(domain[1,2]-domain[2,2]).-1.0

  X1 = chebyshev_polynomial(nodes1,order[1])
  X2 = chebyshev_polynomial(nodes2,order[2])

  X = kron(X2,X1)
  Y = reshape(y,length(nodes1)*length(nodes2),1)

   weights = (X'X)\(X'Y)
   weights = reshape(weights,order[1]+1,order[2]+1)

  return weights

end

function chebyshev_evaluate(w::Array{Float64,2},x::Array{Float64,1},order::Array{Int64,1},domain::Array{Float64,2})

  x1 = 2.0*(x[1]-domain[2,1])/(domain[1,1]-domain[2,1])-1.0
  x2 = 2.0*(x[2]-domain[2,2])/(domain[1,2]-domain[2,2])-1.0

  poly1 = chebyshev_polynomial([x1],order[1])
  poly2 = chebyshev_polynomial([x2],order[2])

  y_hat = 0.0

  for i = 1:order[1]+1
    for j = 1:order[2]+1

      y_hat += w[i,j]*poly1[i]*poly2[j]

    end
  end

  return y_hat

end
