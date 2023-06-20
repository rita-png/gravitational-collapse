using CSV, Tables
using Dierckx


#A = ARGS[1]
#run = ARGS[2]
A=0.0492646484375
run=1
#global dir = "/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/DATA/bisectionsearch"
global dir = ENV["gc_datapath"]#"/home/ritapsantos/data/ritapsantos"

#println(dir)
global bisection = true
global loggrid = true
global meshrefinement = false

include("./ScalarField.jl");

# Parameters

m = 1
res=m;
N=2.0^m*1000.0/2.0#2.0^m*5000.0/2.0#2.0^m*1000.0;#2.0^m*500.0;#N=2.0^m*500.0#2.0^m*100.0;
Xf=1.0;

dx=Xf/N
if loggrid==false
    dt=0.5*round(dx,digits=10)
else
    dt=0.1*round(dx,digits=10)
end
Nt=2.0^m*1000.0/2.0#2.0^m*5000.0/2.0
Tf=Nt*dt;

#### Grid ####

ori=0.0#Float128(0.0)#0.0;
initX1 = nothing
N=int(N)
initX1=range(ori, stop=Xf, step=dx);
initX = range(round(ori-3.0*dx,digits=10), stop=Xf+3.0*dx, step=dx)

L=length(initX);

if loggrid==true
    global originalX=initX
    xtilde=gridfunc(initX1)
    initX1=xtilde
    initX=collect(initX)
    initX[4:L-3]=xtilde
end;

#### Building initial data ####

initm=zeros(L);
initbeta=zeros(L);
initpsi=zeros(L);
initderpsi=zeros(L);

state_array=[initm initbeta initpsi initderpsi initX];

#PSI
r0=0.7#0.3
sigma=0.3
initpsi[4:L-3] = init_gaussian(initX1,r0,sigma,A)

state_array[:,3] = initpsi
state_array = ghost(state_array)

#PSI,X FROM PSI
initderpsi[4:L-3] = init_gaussian_der(initX1,r0,sigma,A)
state_array[:,4] = initderpsi
state_array=ghost(state_array)

####

derpsi_func = Spline1D(initX[4:L-3], state_array[4:L-3,4],  k=4)

# m, beta, psi
y0=[0 0 0]

state_array[4:L-3,1:3] = n_rk4wrapper(RHS,y0,initX[4:L-3],0,derpsi_func,state_array[:,:]);

state_array = ghost(state_array);

run=int(run)
#CSV.write(dir*"/bisectionsearch/run$run/time_step0.csv", Tables.table(state_array), writeheader=false)
###CSV.write("./DATA/bisectionsearch/run$run/time_step0.csv", Tables.table(state_array), writeheader=false)

global files=["m", "beta", "psi", "derpsi"]

print_muninn(files, 0, state_array[:,1:5],res,"w")

time=0
criticality=0.0
explode=0.0
critical_stop=0
evol_stats = [criticality A sigma r0 time explode run]
monitor_ratio = zeros(L)

run=int(run)
if run == 1
    CSV.write(dir*"/bisectionsearch/parameters.csv", Tables.table(evol_stats), writeheader=true, header=["criticality", "A", "sigma", "r0", "time", "explode", "run"])
    #CSV.write("./DATA/bisectionsearch/parameters.csv", Tables.table(evol_stats), writeheader=true, header=["criticality", "A", "sigma", "r0", "time", "explode", "run"])
end

ginit=speed(initX,state_array[:,1],state_array[:,2])

finaltime=1.3
stats,T_interp=timeevolution(state_array,finaltime,run)#timeevolution(state_array,finaltime,dir,run)


CSV.write(dir*"/bisectionsearch/parameters.csv", Tables.table(stats), writeheader=true,header=["criticality", "A", "sigma", "r0", "time", "explode", "run"],append=true);
#CSV.write("./DATA/bisectionsearch/parameters.csv", Tables.table(stats), writeheader=true,header=["criticality", "A", "sigma", "r0", "time", "explode", "run"],append=true);

CSV.write(dir*"/bisectionsearch/timearray.csv", Tables.table(T_interp), writeheader=false);
#CSV.write("./DATA/bisectionsearch/timearray.csv", Tables.table(T_interp), writeheader=false);