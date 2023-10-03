# ---------------------------------------------------------------------------
# Bisection search to look for apparent horizon formation
# ---------------------------------------------------------------------------

using CSV, Tables, DataFrames, Plots, Printf

include("./ScalarField.jl");

# grid 1
"""Agrid=0.35
kgrid=0.7
mgrid=0.5#0.55
fgrid=5"""

#grid 2, newnew
#Agrid=0.42
#kgrid=0.5
#mgrid=0.5
#fgrid=4.5

#grid 3
Agrid=0.44
kgrid=0.47
mgrid=0.5
fgrid=4.0

#grid 4
#Agrid=0.58
#kgrid=0.47
#mgrid=0.51#0.55
#fgrid=1.0

#### CONFIG #### DON'T FORGET TO CHOOSE GRID

global bisection = true
global loggrid = true
global compactified = true
global zeroformat = true


#N=6000.0 #grid 3, tmux 3
#global dir = "/home/ritapsantos/data/ritapsantos/4thuneven"
#global low_bound = 0.1245
#global high_bound = 0.12475#0.1250
#global twod = false

##OR##

#N=10000.0 #grid 2, tmux 2
#global dir = "/home/ritapsantos/data/ritapsantos"
#global low_bound = 0.1245
#global high_bound = 0.12475#0.1250
#global twod = false

##OR##

N=2000.0 #grid 3, tmux 4
global dir = "/home/ritapsantos/data/ritapsantos/new4thuneven"
global low_bound = 0.12474609375000001#0.1245
global high_bound = 0.15#0.12475#0.1250
global twod = false

##OR##
#N=10000 #grid 4, tmux 5
#global dir = "/home/ritapsantos/data/ritapsantos/gridfour"
#global low_bound = 0.1246#0.1245 #these bounds were guessed
#global high_bound = 0.1247#0.1250
#global twod = false
#### CONFIG ####



global run = 1
global runmax = 40

plt_A_crit = Vector{Float64}()
plt_A_non_crit = Vector{Float64}()
plt_x1 = Vector{Int64}()
plt_x2 = Vector{Int64}()


while(run <= runmax)

    A = (low_bound + high_bound) / 2

    println("\n########")
    println("\nBisection search run ##", run, "\nLow bound = ",low_bound,"; High bound = ", high_bound,"; A = ", A," N = ", N,"\n")

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



using Plots

p = scatter(plt_x1, plt_A_non_crit, show=true, xaxis="x", yaxis="Amplitude A", title="Horizon formation",label="Non-critical evolution");
p = scatter!(plt_x2, plt_A_crit,label="Critical evolution");
savefig("bisectionsearch.png")
display(p)

readline()