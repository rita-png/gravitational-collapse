resol=[1]

for ress in resol
    A=0.01
    global ARGS = [A,ress]
    include("./convergencescriptnew.jl");
end
