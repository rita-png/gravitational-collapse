resol=[1,2,3,4]

for ress in resol
    A=0.1
    global ARGS = [A,ress]
    include("./convergencescriptnew.jl");
end
