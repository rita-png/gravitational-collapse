# ---------------------------------------------------------------------------
# Bisection search to look for apparent horizon formation
# ---------------------------------------------------------------------------


#Searching for A* between 0.05 and 0.20

using CSV, Tables, DataFrames, ProgressBars, Plots

global bisection=true
global low_bound = 0.04926157283782959#0.04925#0.04922733163833619
global high_bound = 0.049261573076248164#0.0495#0.04922733211517334
global run = 1
global runmax = 30


plt_A_crit = Vector{Float64}()
plt_A_non_crit = Vector{Float64}()
plt_x1 = Vector{Int64}()
plt_x2 = Vector{Int64}()


while(run <= runmax)

    A = (low_bound + high_bound) / 2

    println("\n########")
    println("\nBisection search run ##", run, "\nLow bound = ",low_bound,"; High bound = ", high_bound,"; A = ", A,"\n")

    global ARGS = [A,run]
    include("./Evolution_ScalarField.jl");
    df = CSV.read("/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/DATA/bisectionsearch/parameters.csv", DataFrame)

    if (df[run+1, :criticality]) == 0.0
        println("\nNon critical!")
        push!(plt_A_non_crit, A)
        push!(plt_x1, run)
    else
        println("\nCritical!")
        push!(plt_A_crit, A)
        push!(plt_x2, run)
    end
    println("\nA = ",df[run+1, :A], " sigma = ", df[run+1, :sigma], " r0 = ", df[run+1, :r0], " Final time = ", df[run+1, :time], " explode = ", df[run+1, :explode])

    #println(df[1, :explode])
    if (df[run+1, :explode]) == 1.0
        println("Found a NaN at time ",df[run+1, :time])
    end
    
    if (df[run+1, :criticality]) == 1.0
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