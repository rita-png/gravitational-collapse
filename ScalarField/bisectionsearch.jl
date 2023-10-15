# ---------------------------------------------------------------------------
# Bisection search to look for apparent horizon formation
# ---------------------------------------------------------------------------

using CSV, Tables, DataFrames, Plots, Printf

include("./ScalarField.jl");

#### CONFIG ####



global N=20000.0
global dir = "/home/ritapsantos/data/ritapsantos/moreres"


####

global bisection = true
global loggrid = false
global compactified = true
global zeroformat = true


#### CONFIG ####



global run = 1
global runmax = 3

plt_A_crit = Vector{Float64}()
plt_A_non_crit = Vector{Float64}()
plt_x1 = Vector{Int64}()
plt_x2 = Vector{Int64}()

As=[0.124650500,0.12465049]
global ii=0
while(run <= runmax)
    global ii=ii+1
    A = As[ii]

    println("\n########")
    println("\nBisection search run ##", run, "; A = ", A," N = ", N,"\n")

    global ARGS = [A,run,N]
    include("./Evolution_ScalarField.jl");
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
    end
    
    if (df[run+1, :Column1]) == 1.0 #:criticality
        global high_bound = A
    else
        global low_bound = A
    end
    global run = run + 1
end