resol=[1,2,3]

for ress in resol
    A=0.05
    global ARGS = [A,ress]
    include("./convergencescriptnew.jl");
end