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

exponents = [-16.25, -15.75, -15.25, -14.75, -14.25, -13.75, -13.25, -12.75, -12.25, -11.75, -11.25, -10.75, -10.25, -9.75, -9.25, -8.75, -8.25, -7.75, -7.25]#collect(-30:0.5:-7)

global run = 1
global runmax = length(exponents)

while(run <= runmax)

    A = exp(exponents[run]) + A_critic

    println("\n########")
    println("\nBisection search run ##", run, "; A = ", A," N = ", N,", exponent = ", exponents[run],"\n")

    global ARGS = [A,run,N]
    include("./Evolution_ScalarField.jl");
    if loggrid==true
        df = CSV.read(dir*"/bisectionsearch/muninnDATA/uneven/parameters.csv", DataFrame)
    else
        df = CSV.read(dir*"/bisectionsearch/muninnDATA/even/parameters.csv", DataFrame)
    end

    println("\nA = ",df[run+1, :Column2], " sigma = ", df[run+1, :Column3], " r0 = ", df[run+1, :Column4], " Final time = ", df[run+1, :Column5], " explode = ", df[run+1, :Column6], " black hole mass = ", df[run+1, :Column7])
    
    
    global run = run + 1
end


