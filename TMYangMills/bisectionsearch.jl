# ---------------------------------------------------------------------------
# Bisection search to look for apparent horizon formation
# ---------------------------------------------------------------------------

using CSV, Tables, DataFrames, Plots, Printf

include("./TMYangMillsTHREE.jl");

#### CONFIG ####

global bisection = true
global loggrid = false
global compactified = true
global zeroformat = true
global twod=false
global source=false
global r0=0.3
global sigma=0.1

####bisec1

#global N=2000.0
#global dir = "/home/ritapsantos/data/ritapsantos/YangMillsbisec1"
#global low_bound=0.088
#global high_bound=0.09


####bisec2

#global N=4000.0
#global dir = "/home/ritapsantos/data/ritapsantos/YangMillsbisec2"
#global high_bound=0.08864099693298341#0.0887
#global low_bound=0.08864097080230714#0.08185

####bisec3

#global N=1000.0
#global dir = "/home/ritapsantos/data/ritapsantos/YangMillsbisec3"
#global low_bound = 0.07
#global high_bound = 0.11

##bisec4

#global N=2000.0
#global dir = "/home/ritapsantos/data/ritapsantos/YangMillsbisec3"
#global low_bound = 0.0885
#global high_bound = 0.08
#

####bisec5

global N=2000.0
global dir = "/home/ritapsantos/data/ritapsantos/YangMillsbisec5"
global high_bound=0.0512109375#0.075
global low_bound=0.0509375#0.04

#### CONFIG ####




global run = 1
global runmax = 60

plt_A_crit = Vector{Float64}()
plt_A_non_crit = Vector{Float64}()
plt_x1 = Vector{Int64}()
plt_x2 = Vector{Int64}()

while(run <= runmax)
    A = (low_bound + high_bound) / 2

    println("\n########")
    println("\nBisection search run ##", run, "\nLow bound = ",low_bound,"; High bound = ", high_bound,"; A = ", A," N = ", N,"\n")

    global ARGS = [A,run,N]

    include("./Evolution_TMTHREEYangMills.jl");
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