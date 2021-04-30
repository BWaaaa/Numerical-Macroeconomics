# Notebook:          III Macro_6 Computational Macro Assignment
# Codes_Main Part I: Solve the labour matching model using PFI
# Author:            2495262
# Date:              Mar 31 2021

using LinearAlgebra
using Statistics
using FastGaussQuadrature
using Plots

include("Tut9_FunctionsII_Newtoniter.jl")
include("Tut9_FunctionsVI_chebyshevandco.jl")

## Set up Parameteres
sigma = 2.00
alpha = 0.6
beta  = 0.99
delta = 0.12
kappa = 0.15
b     = 0.07
m     = 0.4
zeta  = 0.8 # 0.4, 0.5, 0.6, 0.7, 0.8
rho   = 0.95
sd_eps= 0.01

## First we need to solve for the steady state
# We have newton iteration to find it
function lab_match_eqns(x)
    # x is steady state, it includes [n_-1,thet,c; n,thet_+1,c_+1]
    # I can have a dedicated section to derive all the first order conditions
    f = zeros(length(x))
    # Notice the expectation is non-relevant
    # This is because at steady state, we assume no shock (expz=1)
    # Therefore the expectation is void, i.e., the averaging
    # just gives the actual value
    f[1] = (1-delta)*x[1] + m*(1-(1-delta)*x[1])*abs(x[2])^(1-alpha)-x[4]
    # this equation indicates the new unemployement is a process
    # of those lose job and those not finding jobs
    # therefore, the employment n is one minus, plus job formation
    # Abs is included to avoid negative exponential issues (e.g. (-0.2)^0.4)
    f[2] = x[3] + kappa*(1-(1-delta)*x[1])*abs(x[2]) - exp(0)*x[4]
    # the LHS is output, and it is consumed or burnt because of
    # job formation---frivalous frictional job finding
    # We have no capital in this model! great simplication.
    # Therefore, we do not have a part for savings
    f[3] = (1-zeta)*(exp(0)-b) + beta*(1-delta)* x[6]^(-sigma) *
    (1-zeta*m*x[5]^(1-alpha))/x[3]^(-sigma)*kappa/m*x[5]^alpha -
    kappa/m*abs(x[2])^alpha
    # This is the relatively complicated wage determination process
    # which is determined by Nash Bargaining parameter zeta
    # In the report, I will present the derivation of this equation
    # and discuss the meaning of Nash Bargaining with zeta
    # It is made up of two parts: unemployment benefit,
    # plus discounted continuation value
    f[4] = x[4] - x[1]
    f[5] = x[5] - x[2]
    f[6] = x[6] - x[3]
    return f
end
tol      = 1e-8
maxiter1 = 100
initval  = [0.85, 10.0, 1.5, 0.85, 10.0, 1.5]
(ss,f_ss,iters)=newton_iter(lab_match_eqns, initval, tol, maxiter1)

## Next, we proceed to solving the PFI problem
# we, again, have two states, technology and labour
num_nodes_tech = 9
order_tech     = 5
tech_low       = sqrt(sd_eps^2/(1-rho^2))*(-4.0)
tech_high      = sqrt(sd_eps^2/(1-rho^2))*4.0
domain_tech    = [tech_high,tech_low]
nodes_tech     = chebyshev_nodes(num_nodes_tech, domain_tech)

num_nodes_lab  = 15
order_lab      = 6
lab_low        = ss[1]*0.5 # or 0.8, adjust is essential
lab_high       = 1.0 # the mploymet, by assumption of unit mass and its interpretation as a rate, <1
domain_lab     = [lab_high, lab_low]
nodes_lab      = chebyshev_nodes(num_nodes_lab, domain_lab)

order = [order_tech, order_lab]
domain = [domain_tech domain_lab]
eps_nodes, eps_weights = gausshermite(11)

## Third, we initialize the policy function
# notice, the policy function is a triplet: (n_t, c_t, thet_t)
# We are solving a three-equation three unknown system
# Intake? States = [tech, n_t-1]
labour      = ones(num_nodes_tech, num_nodes_lab)*ss[1]
tightness   = ones(num_nodes_tech, num_nodes_lab)*ss[2]
consumption = ones(num_nodes_tech, num_nodes_lab)*ss[3]

# weights, although seems like coefficients, are the function we appx
# We have two function so approximate, becase for labour, we do not
# need to calculate future avergaing, and thus cleanly separated
consumption_weights = chebyshev_weights(consumption, nodes_tech,nodes_lab,order,domain)
tightness_weights   = chebyshev_weights(tightness, nodes_tech,nodes_lab,order,domain)

init = [labour[1,1], tightness[1,1], consumption[1,1]]

## Fourth, we show What we do for each state, by doing the first iteration's first state
state = [nodes_tech[1],nodes_lab[1]]
function stoch_labmatch_model(x)
    # This is to calculate the value of difference and the aim is to converge 0
    # I am saying: given the policy, how far am I away from fix point
    # So: x = [n, thet, c]
    f = zeros(length(x))
    f[1] = (1-delta)*state[2] + m*(1-(1-delta)*state[2])*abs(x[2])^(1-alpha) - x[1]
    # This is indeed the process solver: given states, and i have
    # thet_t coming from x as initial value, i calculate n and report distance
    # Abs is included to avoid negative exponential issues (e.g. (-0.2)^0.4)
    f[2] = x[3] + kappa*(1-(1-delta)*state[2])*x[2] - exp(state[1])*x[1]
    future_tech_nodes = rho*state[1] .+ sqrt(2)*sd_eps*eps_nodes
    # Given future_tech_nodes, we have state for next period
    # state_+1 = [future_tech, x[1]], in which x[1]=n
    # Now, we can construct the continue_value
    continue_val = 0.0
    state_new    = zeros(2)
    # This is to plug in for finding c+1 and thet_+1
    # Notice that they are both functions, and we just need state to input
    # The function, of course, are the weights
    for i = 1:11
        state_new    .= [future_tech_nodes[i], x[1]]
        thet_future   = chebyshev_evaluate(tightness_weights, state_new, order, domain)
        continue_val += pi^(-1/2)*eps_weights[i] * 1/abs(x[3])^(-sigma) * kappa/m *
        abs(chebyshev_evaluate(consumption_weights, state_new, order, domain))^(-sigma) *
        (1-zeta*m*abs(thet_future)^(1-alpha))*abs(thet_future)^alpha
        # this is using the chebyshev weights for con and tight to
        # calculate the integration
    end
    f[3] = (1-zeta)*(exp(state[1])-b) + beta*(1-delta)*continue_val - kappa/m*abs(x[2])^alpha
    # For this one, we need to calculate continuation value
    # Notice in this "euler", we divided c_t from LHS to RHS, and thus
    # there is a pert x[3] not stochastic in the integration
    return f
end
first_state_first_iter = stoch_labmatch_model(init)
tol                    = 1e-6
maxiters               = 500
(soln_1, f_soln, iters)= newton_iter(stoch_labmatch_model, init, tol, maxiters)

## Fiifth, we start iterating on the three-equatio npolicy function structure
new_lab     = similar(labour)
new_tight   = similar(tightness)
new_con     = similar(consumption)
while true
    consumption_weights .= chebyshev_weights(consumption, nodes_tech, nodes_lab, order, domain)
    tightness_weights   .= chebyshev_weights(tightness, nodes_tech, nodes_lab, order, domain)
    # Form the functions, based on our last iteratio, this step
    # is closing our target last round: respect the function al form
    # Otherwise, declare global
    for i = 1:num_nodes_tech
        for j = 1:num_nodes_lab
            # Iterating over all grid points
            state .= [nodes_tech[i], nodes_lab[j]]
            init  .= [labour[i,j], tightness[i,j], consumption[i,j]]
            (soln,f_soln,iters) = newton_iter(stoch_labmatch_model, init, tol, maxiters)
            # Notice I print the number of iterations
            new_lab[i,j]   = soln[1]
            new_tight[i,j] = max(soln[2],1e-3)
            # to ensure that varibles are always positive strictly
            new_con[i,j]   = max(soln[3],1e-3)
        end
    end
    len = maximum(abs, [(new_lab-labour) (new_tight-tightness) (new_con-consumption)])
    println("\nlen=",len)
    if len<tol
        break
    end
    labour      .= new_lab
    tightness   .= new_tight
    consumption .= new_con
end

## Finally, let's write the results into an output file
suffix   = Int(floor(mod(time(),100)*10))
filename = "PFIResults_zeta0"*string(Int(zeta*10))*"_time"*string(suffix)*".txt"
path     = @__DIR__
# or yourdirectory
cd(path)
open(filename,"w") do io
    write(io,"zeta is ")
    println(io, zeta)
    write(io, "\nthe steady state is (n,thet,c)=")
    println(io, ss[1:3])
    write(io,"\nGuide: each roll reads: given same tech, what is the policy at differnt lab_-1\n\n")
    write(io, "Labour(tech,lab_-1)=size")
    println(io,size(labour))
    for li in 1:size(labour)[1]
        println(io, labour[li,:])
    end
    write(io, "\nTightness(tech,lab_-1)=size")
    println(io,size(tightness))
    for li in 1:size(tightness)[1]
        println(io, tightness[li,:])
    end
    #write(io, "Tightness")
    write(io, "\nConsumption(tech,lab_-1)=size")
    println(io,size(consumption))
    for li in 1:size(consumption)[1]
        println(io, consumption[li,:])
    end
end

## As a preperation for next analysis, let us plot the policy functions
figpath = path*"/figures/" #"\\figures\\"
# or your directory
if !ispath(figpath)
    mkdir(figpath)
end
function plot_foo()
    mid = Int(floor((num_nodes_tech+1)/2))
    f1 = plot(nodes_lab,[labour[1,:], labour[mid,:], labour[end,:]], label=["low tech" "mid tech" "high tech"], legend = :topleft)
    #savefig(figpath*"Tut8_policyn_zeta0"*string(Int(zeta*10))*"_time"*string(suffix))
    f2 = plot(nodes_lab,[tightness[1,:], tightness[mid,:], tightness[end,:]], label=["low tech" "mid tech" "high tech"], legend = :topleft)
    #savefig(figpath*"Tut8_policythet_zeta0"*string(Int(zeta*10))*"_time"*string(suffix))
    f3 = plot(nodes_lab,[consumption[1,:], consumption[mid,:], consumption[end,:]], label=["low tech" "mid tech" "high tech"], legend = :topleft)
    #savefig(figpath*"Tut8_policyc_zeta0"*string(Int(zeta*10))*"_time"*string(suffix))
    f4 = plot(f1, f2, f3, layout=(1,3))
    return f4
end
f5 = plot_foo()
savefig(figpath*"/PolicyFunctions_zeta0"*string(Int(zeta*10))*"_time"*string(suffix))
