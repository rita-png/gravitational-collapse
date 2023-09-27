# ---------------------------------------------------------------------------
# Bisection search to look for apparent horizon formation
# ---------------------------------------------------------------------------

using CSV, Tables, DataFrames, Plots, Printf

include("./ScalarField.jl");

#### CONFIG ####

####started 21/09####
global N=2000.0
global dir = "/home/ritapsantos/data/ritapsantos/massscaling"



####

global bisection = true
global loggrid = false
global compactified = true
global zeroformat = true
global twod = false

#### CONFIG ####


A_critic = 0.12465049029170887

exponents = collect(-30:0.5:-15)

global run = 1
global runmax = length(exponents)

while(run <= runmax)

    A = exp(exponents[run]) + A_critic

    println("\n########")
    println("\nBisection search run ##", run, "\nLow bound = ",low_bound,"; High bound = ", high_bound,"; A = ", A," N = ", N,"\n")

    global ARGS = [A,run,N]
    """include("./Evolution_ScalarField.jl");
    if loggrid==true
        df = CSV.read(dir*"/bisectionsearch/muninnDATA/uneven/parameters.csv", DataFrame)
    else
        df = CSV.read(dir*"/bisectionsearch/muninnDATA/even/parameters.csv", DataFrame)
    end

    if (df[run+1, :Column1]) == 0.0 #:criticality
        println("\nNon critical!")
        push!(plt_A_non_crit, A)
        push!(plt_x1, run)
    else
        println("\nCritical!")
        push!(plt_A_crit, A)
        push!(plt_x2, run)
    end
    println("\nA = ",df[run+1, :Column2], " sigma = ", df[run+1, :Column3], " r0 = ", df[run+1, :Column4], " Final time = ", df[run+1, :Column5], " explode = ", df[run+1, :Column6])

    #println(df[1, :explode])
    if (df[run+1, :Column6]) == 1.0 #:explode
        println("Found a NaN")
    end"""
    
    
    global run = run + 1
end


