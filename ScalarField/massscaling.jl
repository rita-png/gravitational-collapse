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


#A_critic = 0.12465049029170887 #oolon
A_critic = 0.12465049029170887 #baltasar


exponents = collect(-30:0.5:-7)

global run = 1
global runmax = length(exponents)

while(run <= runmax)

    A = exp(exponents[run]) + A_critic

    println("\n########")
    println("\nBisection search run ##", run, "; A = ", A," N = ", N,", exponent = ", exponents[run],"\n")

    global ARGS = [A,run,N,sigma]
    include("./Evolution_ScalarField_Diff_Families.jl");
    if loggrid==true
        df = CSV.read(dir*"/bisectionsearch/muninnDATA/uneven/parameters.csv", DataFrame)
    else
        df = CSV.read(dir*"/bisectionsearch/muninnDATA/even/parameters.csv", DataFrame)
    end

    println("\nA = ",df[run+1, :Column2], " sigma = ", df[run+1, :Column3], " r0 = ", df[run+1, :Column4], " Final time = ", df[run+1, :Column5], " explode = ", df[run+1, :Column6], " black hole mass = ", df[run+1, :Column7])
    
    
    global run = run + 1
end


