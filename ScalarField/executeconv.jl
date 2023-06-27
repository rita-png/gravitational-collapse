resol=[1,2,3]

for ress in resol
    A=0.01
    N=100
    global ARGS = [A,ress,N]
    include("./convergencescriptnew.jl");
end