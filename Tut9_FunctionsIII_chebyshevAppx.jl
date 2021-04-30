# Notebook:            Lecture 5 Function Approximation
# Codes_Functions III: Chebyshev Polynomials
# Author:              2495262
# Date:                Mar 13 2021

function polynomial_approx(y_points,x_points,point,n)

  x_lagrange = zeros(length(x_points),n+1)
  for i = 1:(n+1)

    x_lagrange[:,i] .= x_points.^(i-1)

  end

  alpha = (x_lagrange'x_lagrange)\(x_lagrange'y_points)

  y_hat = 0.0
  for i in 1:n+1
    y_hat += alpha[i]*point^(i-1)
  end

  return y_hat

end


function chebyshev_nodes(N)

  x = zeros(N)
  for i = 1:N

    x[i] = - cos((2*i-1)*pi/(2*N))

  end

  return x

end


function chebyshev_nodes(N::Integer,domain::Array{T,1}) where {T <: Real}

  #domain = [upper,lower]

  x = chebyshev_nodes(N) # between [-1,1]

  # when x[i] = -1 then x_scaled[i] = lower
  # when x[i] = 1 then x_scaled[i] = upper

  x_scaled = (domain[1]+domain[2])/2 .+ ((domain[1]-domain[2])/2)*x

  return x_scaled

end


function chebyshev_polynomial(x_points,n)

  poly_matrix = zeros(length(x_points),n+1)

  for i = 1:(n+1)

    if i == 1
      poly_matrix[:,i] .= ones(length(x_points))
    elseif i == 2
      poly_matrix[:,i] .= x_points
    else
      poly_matrix[:,i] .= 2*x_points.*poly_matrix[:,i-1] - poly_matrix[:,i-2]
    end

  end

  return poly_matrix

end
