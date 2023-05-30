using Interpolations
using ProgressMeter
#using Interpolations
using PyCall
using CSV, Tables
using Dierckx
using Printf

"""println("#######3")
println("argument received is", ARGS[1])
println("argument type received is", typeof(ARGS[1]))
println("#######3")"""



A = ARGS[1]
run = ARGS[2]

include("./ScalarField.jl");


# Parameters


m = 1
res=m;
N=2.0^m*7500.0/2.0#2.0^m*2000.0/2.0
Xf=1.0;

dx=Xf/N;
dt=0.5*round(dx,digits=10);
Nt=2.0^m*7500.0/2.0
Tf=Nt*dt; #final time


global dir = "/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/DATA/bisectionsearch"

ori=0.0
initX1 = nothing
N=int(N)
initX1=range(ori, stop=Xf, step=dx);
#initX1=create_range(ori,Xf,dx,N)
initX = range(round(ori-3.0*dx,digits=10), stop=Xf+3.0*dx, step=dx)
#initX=create_range(ori-3.0*dx,Xf+3.0*dx,dx,N+6)

L=length(initX);

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

derpsi_func = Spline1D(initX[4:L-3], initderpsi[4:L-3],  k=4)

y0=[0.0 0.0 0.0]
#state_array[4:L-3,1:3] = twod_n_rk4wrapper(RHS,y0,initX[4:L-3],0,derpsi_func,state_array[:,:],coef);
state_array[4:L-3,1:3] = twod_n_rk4wrapper(RHS,y0,initX[4:L-3],0,derpsi_func,state_array[:,:]);

state_array = ghost(state_array);

run=int(run)
CSV.write(dir*"/run$run/time_step0.csv", Tables.table(state_array), writeheader=false)

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
    CSV.write(dir*"/parameters.csv", Tables.table(evol_stats), writeheader=true, header=["criticality", "A", "sigma", "r0", "time", "explode", "run"])
end

ginit=speed(initX,state_array[:,1],state_array[:,2],dx)
#println(update_dt(initX,state_array[:,1],state_array[:,2],dx,ginit)/dt)

finaltime=1.3
stats,T_interp=timeevolution(state_array,finaltime,dir,run)


CSV.write("/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/DATA/bisectionsearch/parameters.csv", Tables.table(stats), writeheader=true,header=["criticality", "A", "sigma", "r0", "time", "explode", "run"],append=true);

CSV.write("/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/DATA/bisectionsearch/timearray.csv", Tables.table(T_interp), writeheader=false,append=true);