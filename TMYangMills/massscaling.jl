# ---------------------------------------------------------------------------
# Bisection search to look for apparent horizon formation
# ---------------------------------------------------------------------------

using CSV, Tables, DataFrames, Plots, Printf

include("./TMYangMillsTHREE.jl");

#### CONFIG ####

global N=4000.0
global dir = "/home/ritapsantos/data/ritapsantos/YMmassscaling"


####

global bisection = true
global loggrid = false
global compactified = true
global zeroformat = true
global twod = false
global source=false
global r0=0.3
global sigma=0.1

#### CONFIG ####


A_critic = 0.0886409955039620

exponents = collect(-25:0.5:-15.5)#collect(-15:0.5:-1)#collect(-25:0.5:-19)#collect(-19.5:0.5:-15.5)#collect(-15:0.5:-4.5)#[-18.5,-18,-17.5]#[-9.5, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5,  0.0]#[-26,-25.5,-25,-24.5,-24,-23.5,-23,-22.5,-22,21.5,-21,-20.5,-20,-19.5,-19,-18.5,-18,-17.5,17,-16.5]#collect(-17:0.5:-10)#collect(-30:0.25:-7)
global run = 1
global runmax = length(exponents)

while(run <= runmax)

    A = -exp(exponents[run]) + A_critic

    println("\n########")
    println("\nBisection search run ##", run, "; A = ", A," N = ", N,", exponent = ", exponents[run],"\n")

    global ARGS = [A,run,N,sigma,r0]
    include("./Evolution_YangMills_Diff_Families.jl");
    if loggrid==true
        df = CSV.read(dir*"/bisectionsearch/muninnDATA/uneven/parameters.csv", DataFrame)
    else
        df = CSV.read(dir*"/bisectionsearch/muninnDATA/even/parameters.csv", DataFrame)
    end

    println("\nA = ",df[run+1, :Column2], " sigma = ", df[run+1, :Column3], " r0 = ", df[run+1, :Column4], " Final time = ", df[run+1, :Column5], " explode = ", df[run+1, :Column6], " black hole mass = ", df[run+1, :Column7])
    
    
    global run = run + 1
end


