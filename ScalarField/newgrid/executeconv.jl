resol=[4,5]

for ress in resol
    A=0.1#0.01
    global ARGS = [A,ress]
    include("./convergencescriptnew.jl");
end