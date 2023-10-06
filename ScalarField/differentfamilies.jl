# ---------------------------------------------------------------------------
# Bisection search to look for apparent horizon formation
# ---------------------------------------------------------------------------

using CSV, Tables, DataFrames, Plots, Printf

include("./ScalarField.jl");

global bisection = true
global loggrid = true
global compactified = true
global zeroformat = false



N=2000.0
global subdir = "/home/ritapsantos/data/ritapsantos/differentfamilies"
global low_bound = 0.08
global high_bound = 0.170
global twod = false


#find A* for different sigmas than 0.3 - doing this rn
#or find sigma* for an A*?

global run = 1
global runmax = 30
sigmas = [0.1 0.2 0.4 0.5 0.6]

global kk=0
for sigma in sigmas
    global kk=kk+1
    global dir = subdir*string(kk)
    println("dir is ", dir)
    while(run <= runmax)

        A = (low_bound + high_bound) / 2

        println("\n########")
        println("\nBisection search for sigma = ", sigma,"; run ##", run, "\nLow bound = ",low_bound,"; High bound = ", high_bound,"; A = ", A," N = ", N,"\n")

        global ARGS = [A,run,N,sigma]
        include("./Evolution_ScalarField_Diff_Families.jl");
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

end