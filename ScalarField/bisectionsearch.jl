# ---------------------------------------------------------------------------
# Bisection search to look for apparent horizon formation
# ---------------------------------------------------------------------------

using CSV, Tables, DataFrames, Plots, Printf

include("./ScalarField.jl");

#### CONFIG ####

"""global N=6000.0
global dir = "/home/ritapsantos/data/ritapsantos/4theven"

global low_bound = 0.12465049985278165#0.12465049985272345#0.12465049985051155#0.12465049743652343#0.12465045166015624
global high_bound = 0.12465049985275255#0.12465049985283985#0.12465049985796213#0.1246505012512207#0.1246505126953125

"""
##OR##
####ongoing####
#global N=20000.0
#global dir = "/home/ritapsantos/data/ritapsantos"
#
##0.124648615151412<A*<0.125#
#global low_bound = 0.124648615151412#0.124#0.12465050144332973#0.12465050144332965#0.12465049985051155
#global high_bound = 0.125#0.12465052097457974#0.12466050144332981#0.12465050251698608#0.125

##OR##
####started 21/09####
global N=2000.0
global dir = "/home/ritapsantos/data/ritapsantos/massscaling"

global low_bound = 0.124
global high_bound = 0.125
####

global bisection = true
global loggrid = false
global compactified = true
global zeroformat = true
global twod = false


#global dir = "/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/DATA"

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