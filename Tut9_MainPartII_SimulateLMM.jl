# Notebook:           III Macro_6 Computational Macro Assignment
# Codes_Main Part II: Simulate the Labour Macthing Model, focusing on consumption
# Author:             2495262
# Date:               Apr 3 2021

using LinearAlgebra
using Random
using Statistics
using Plots

include("Tut9_FunctionsVI_chebyshevandco.jl")

## Set up Parameteres and read the input file
# and read the file
path = @__DIR__
# or your path
cd(path)
filename = "PFIResults_zeta08_time802.txt"
# please input the filename you just created in Codes_Main Part I
open(filename, "r") do io
    line  = readline(io)
    line1 = match(r"[0-9].[0-9]",line)
    global zeta = parse(Float64,line1.match)
    # global helps to store value in global
end
sigma = 2.00
alpha = 0.6
beta  = 0.99
delta = 0.12
kappa = 0.15
b     = 0.07
m     = 0.4
# zeta as read from the input file
rho   = 0.95
sd_eps= 0.01

## Read the steady state
function read_steadystate(filenm)
    store = zeros(3)
    open(filenm, "r") do io
        line1  = readline(io)
        line2  = readline(io)
        line3  = readline(io)
        seec   = SubString.(line3, findall(r"[0-9]+.[0-9]+",line3))
        store .= parse.(Float64,seec)
        # global helps to store value in global
    end
    return store
end
ss    = read_steadystate(filename)
n0    = ss[1]
tech0 = 0
# we want the economy to start from the steady state
# so similarly tech0 = 0

## and read the policy function
function read_policy(filenm)
    open(filenm, "r") do io
        readline(io), readline(io), readline(io)
        readline(io), readline(io), readline(io)
        line7  = readline(io)
        seec   = SubString.(line7, findall(r"[1-9]+",line7))
        global dim = (parse(Int8,seec[2]), parse(Int8, seec[3]))
        # global helps to store value in global
    end
    store = zeros(3*dim[1], dim[2])
    idx   = [0]
    for line in eachline(filenm)
        if  match(r"\[", line) != nothing
            println("good line", line)
            idx[1]  = idx[1] + 1
            if idx[1] >= 2
                rule         = r"[0-9]+.[0-9]+"
                seec         = SubString.(line, findall(rule,line))
                store[idx[1]-1,:] = parse.(Float64, seec)
            end
        end
    end
    labour      = store[1:dim[1], :]
    tightness   = store[dim[1]+1:2*dim[1],:]
    consumption = store[2*dim[1]+1:3*dim[1],:]
    return labour, tightness, consumption, dim
end
order_tech     = 5
tech_low       = sqrt(sd_eps^2/(1-rho^2))*(-4.0)
tech_high      = sqrt(sd_eps^2/(1-rho^2))*4.0
domain_tech    = [tech_high,tech_low]
order_lab      = 6
lab_low        = ss[1]*0.5 # 0.8, adjust is essential
lab_high       = 1.0 # the mploymet, by assumption of unit mass and its interpretation as a rate, <1
domain_lab     = [lab_high, lab_low]
order = [order_tech, order_lab]
domain = [domain_tech domain_lab]
function approximate_policy(lab, tight, con, dim)
    num_nodes_tech = dim[1]
    nodes_tech     = chebyshev_nodes(num_nodes_tech, domain_tech)
    num_nodes_lab  = dim[2]
    nodes_lab      = chebyshev_nodes(num_nodes_lab, domain_lab)

    labour_weights      = chebyshev_weights(lab, nodes_tech,nodes_lab,order,domain)
    tightness_weights   = chebyshev_weights(tight, nodes_tech,nodes_lab,order,domain)
    consumption_weights = chebyshev_weights(con, nodes_tech,nodes_lab,order,domain)
    return labour_weights, tightness_weights, consumption_weights
end
lab_policy, tight_policy, con_policy, dim = read_policy(filename)
labour_weights, tightness_weights, consumption_weights = approximate_policy(lab_policy, tight_policy, con_policy, dim)

## And with all the input, we complete the economy profile with tech+5+thet varaibles
# Economy is described by [tech, employment, tightness, consumption, output, unemployment, vacancies]
function ecoonomy_vars(tech,n0, tt)
    # tech is the simulation of technology process
    # We also need an initial, steady state level employment n0
    n, thet, c, Y, u, v = zeros(tt), zeros(tt), zeros(tt), zeros(tt), zeros(tt), zeros(tt)
    state   = [tech[1], n0]
    n[1]    = chebyshev_evaluate(labour_weights, state, order, domain)
    thet[1] = chebyshev_evaluate(tightness_weights, state, order, domain)
    c[1]    = chebyshev_evaluate(consumption_weights, state, order, domain)
    Y[1]    = exp(state[1])*n[1]
    u[1]    = 1-(1-delta)*state[2]
    v[1]    = thet[1]*u[1]
    for t in 2:tt
        state .= [tech[t],n[t-1]]
        # println("employment at time", t, state)
        n[t]    = chebyshev_evaluate(labour_weights, state, order, domain)
        thet[t] = chebyshev_evaluate(tightness_weights, state, order, domain)
        c[t]    = chebyshev_evaluate(consumption_weights, state, order, domain)
        Y[t]    = exp(state[1])*n[t]
        u[t]    = 1-(1-delta)*state[2]
        v[t]    = thet[t]*u[t]
    end
    return n, thet, c, Y, u, v
end
# To calculate welfare, we use utility function
function welfare(x)
    # x = [c1,c2,c3,..., ct] is a series of consumption decision
    per     = length(x)
    welfare = 0.0
    for t in 1:per
        welfare += beta^(t-1)*(x[t]^(1-sigma)-1)/(1-sigma)
    end
    return welfare
end

## Now implement the above simulation for one, and then 100_000 Times
ts_per     = 50
TimeSeries = 100_000
# for each simulation, we have one welfare result
Wfare = zeros(TimeSeries)
# for each simluation, we have a time series of [n,thet,c,Y,u,v], we are interested in the mean
mean_employment, mean_tightness, mean_consumption, mean_output, mean_unemployment, mean_vacancy = zeros(TimeSeries), zeros(TimeSeries), zeros(TimeSeries), zeros(TimeSeries), zeros(TimeSeries), zeros(TimeSeries)

Random.seed!(5086_2495262_20210331)
eps   = randn(TimeSeries, ts_per)*sd_eps
for N in 1:TimeSeries
    tech    = zeros(ts_per)
    tech[1] = rho*tech0 + eps[N,1]
    for t in 2:ts_per
        tech[t] = rho*tech[t-1] + eps[N,t]
    end
    employment, tightness, consumption, Y_output, unemployment, vacancy = ecoonomy_vars(tech, n0, ts_per)
    mean_employment[N], mean_tightness[N], mean_consumption[N], mean_output[N], mean_unemployment[N], mean_vacancy[N] = mean(employment), mean(tightness), mean(consumption), mean(Y_output), mean(unemployment), mean(vacancy)
    #println(tech[1], employment[1], tightness[1], consumption[1], Y_output[1], unemployment[1], vacancy[1])
    Wfare[N] = welfare(consumption)
    if mod(N,1_000) == 0
        println("Simulation", N," welfare=", Wfare[N])
    end
end
mean_w = mean(Wfare)

## Create figures of mean_employment, mean_consumption and welfare
figpath = path*"/figures/" #"\\figures\\"
# or your directory
if !ispath(figpath)
    mkdir(figpath)
end
function plot_ncW_foo()
    h1 = histogram(mean_employment, bin=100, label="mean employment, \nwhere steady state is\n"*string(ss[1]))
    # the height is not 1, it is because normalize is to ensure the
    # area under the curve is one, not to say each rectangle has heigh 1
    h2 = histogram(mean_consumption, bin=100, label="mean consumption, \nwhere steady state is\n"*string(ss[3]))
    h3 = histogram(Wfare, bin=100, label="welfare" )
    # the value of welfare depends on the choice of ts_oer
    h4 = plot(h1,h2,h3,layout=(1,3))
    return h4
end
h5 = plot_ncW_foo()
savefig(figpath*"MeanKeyVars_Welfare_dis_zeta0"*string(Int(zeta*10))*"ts"*string(ts_per))

# Create figuresfor other economic varaibles on demand
function plot_allvar_foo()
    h1 = histogram(mean_employment, bin=100, label="mean employment, \nwhere steady state is\n"*string(ss[1]))
    h2 = histogram(mean_consumption, bin=100, label="mean consumption, \nwhere steady state is\n"*string(ss[3]))
    h3 = histogram(mean_tightness, bin=100, label="mean tightness, \nwhere steady state is\n"*string(ss[2]))
    h4 = histogram(mean_output, bin=100, label="mean output")
    h5 = histogram(mean_unemployment, bin=100, label="mean unemployment")
    # notice delta = 0.12, so mean_unemployment + mean_employment = 1.12
    h6 = histogram(mean_vacancy, bin=100, label="mean vacancy")
    h7 = plot(h1,h2,h3,h4,h5,h6,layout=(2,3))
    return h7
end
h8 = plot_allvar_foo()
savefig(figpath*"MeanAllVar_dis_zeta0"*string(Int(zeta*10))*"ts"*string(ts_per))

## For one simulation of time series, we show the exact process below
# And use it to create timespan figures
Random.seed!(5086_2495262_20210331)
eps     = randn(ts_per)*sd_eps
tech    = zeros(ts_per)
tech[1] = rho*tech0 + eps[1]
for t in 2:ts_per
    tech[t] = rho*tech[t-1] + eps[t]
end
employment, tightness, consumption, Y_output, unemployment, vacancy = ecoonomy_vars(tech, n0, ts_per)

function plot_ts_foo()
    p1 = plot(1:ts_per, tech, label="technology process")
    p2 = plot(1:ts_per,[employment, ones(ts_per)*ss[1]], label=["employment" "steady state value"])
    p3 = plot(1:ts_per,[tightness, ones(ts_per)*ss[2]], label=["tightness" "steady state value"])
    p4 = plot(1:ts_per,[consumption, ones(ts_per)*ss[3]], label=["consumption" "steady state value"])
    p5 = plot(p1,p2,p3,p4,layout=(2,2))
    return p5
end
p6 = plot_ts_foo()
savefig(figpath*"KeyVar_timespan_zeta0"*string(Int(zeta*10)))
