# ---------------------------------------------------------------------------
# Bisection search to look for apparent horizon formation
# ---------------------------------------------------------------------------

using CSV, Tables, DataFrames, Plots, Printf

include("./ScalarField.jl");




#### CONFIG ####

#2nd
global low_bound = 0.09#0.10#0.1244#0.140#0.14#0.15109374999#0.15109375
global high_bound = 0.12440009765624999#0.1246#0.141#0.16#0.15109375#0.15109375000698494#0.15121093749999998#0.15132812499999998#0.1515625

## 4th
#global low_bound = 0.151529541015625#0.15151855468749997
#global high_bound = 0.1515350341796875#0.15154052734375#0.1515625

global bisection = true
global loggrid = true
global compactified = true
global zeroformat = true
global twod = true

#global dir = "/home/ritapsantos/data/ritapsantos/4theven"

if twod==false
    global dir = "/home/ritapsantos/data/ritapsantos"
else
    global dir = "/home/ritapsantos/data/ritapsantos/2nduneven"
end
#global dir = "/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/DATA"

if twod==false
    global N=800.0
else
    global N=800.0 #N=1000.0
end

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