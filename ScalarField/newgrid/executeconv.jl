resol=[1]

for ress in resol
    A=0.01#0.1
    global ARGS = [A,ress]
    include("./convergencescriptnew.jl");
end