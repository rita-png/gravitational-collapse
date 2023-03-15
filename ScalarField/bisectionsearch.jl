# ---------------------------------------------------------------------------
# Bisection search to look for apparent horizon formation
# ---------------------------------------------------------------------------


#Searching for A* between 0.05 and 0.20

using CSV, Tables, DataFrames, ProgressBars, Plots


global low_bound = 0.05
global high_bound = 0.20
global run = 1
global runmax = 5


plt_A_crit = Vector{Float64}()
plt_A_non_crit = Vector{Float64}()


while(run <= runmax)

    A = (low_bound + high_bound) / 2

    println("\n########")
    println("\nBisection search run ##", run, "\nLow bound = ",low_bound,"; High bound = ", high_bound,"; A = ", A,"\n")

    global ARGS = [A]
    include("./Evolution_ScalarField.jl");
    df = CSV.read("./DATA/bisectionsearch/parameters.csv", DataFrame)

    if (df[1, :criticality]) == 0.0
        println("\nNon critical!")
        push!(plt_A_non_crit, A)
        push!(plt_x1, A)
    else
        println("\nCritical!")
        push!(plt_A_crit, A)
        push!(plt_x2, A)
    end
    println("\nA = ",df[1, :A], " sigma = ", df[1, :sigma], " r0 = ", df[1, :r0], " Final timestep = ", df[1, :timestep], " explode = ", df[1, :explode])

    #println(df[1, :explode])
    if (df[1, :explode]) == 1.0
        println("Found a NaN at timestep ",df[1, :timestep])
    end
    
    if (df[1, :criticality]) == 1.0
        global high_bound = A
    else
        global low_bound = A
    end
    global run = run + 1
end



using Plots


plt_x1= [1,2,3]
plt_x2= [4,5]

plt_y1=[1, 4, 6]
plt_y2=[1, 1]

p = scatter(plt_x1, plt_y1, show=true, xaxis="x", yaxis="Amplitude A", title="Horizon formation",label="Critical evolution");
p = scatter!(plt_x2, plt_y2,label="Non-critical evolution");
display(p)

readline()