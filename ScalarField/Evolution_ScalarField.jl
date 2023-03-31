using Interpolations
using ProgressMeter
using Interpolations
using PyCall
using CSV, Tables

"""println("#######3")
println("argument received is", ARGS[1])
println("argument type received is", typeof(ARGS[1]))
println("#######3")"""



A = ARGS[1]
run = ARGS[2]

include("./ScalarField.jl");
scipy = pyimport("scipy")
scipyinterpolate = pyimport("scipy.interpolate")


# Parameters
"""m = 1
res=m;
N=100.0#2.0^m*1000.0;#2.0^m*500.0;#N=2.0^m*500.0#2.0^m*100.0;
Xf=1.0;

dx=Xf/(2.0^m*500.0)
dt=round(dx,digits=10);
Nt=100.0#100.0*2^m*10
Tf=1#Nt*dt; #final time"""

m = 1
res=m;
N=2.0^m*750.0/2.0#2.0^m*1000.0;#2.0^m*500.0;#N=2.0^m*500.0#2.0^m*100.0;
Xf=1.0;

dx=Xf/N;
dt=round(dx,digits=10);
Nt=2^m*750.0/2.0#100.0*2^m*10
Tf=1#Nt*dt; #final time
#print(Tf)

# Setting RESOLUTION

global dir = "/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/DATA/bisectionsearch"

ori=0.0;
initX1 = nothing
initX1=range(ori, stop=Xf, step=dx);

initX = range(round(ori-3.0*dx,digits=10), stop=Xf+3.0*dx, step=dx)

L=length(initX);
T=range(0,stop=Tf,step=dt)


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
"""for i in 4:L-3
    initderpsi[i]=Der(state_array[:,:],i,3,initX)
end

state_array[:,4] = initderpsi
state_array = ghost(state_array)"""

#new
initderpsi[4:L-3] = init_gaussian_der(initX1,r0,sigma,A)
state_array[:,4] = initderpsi
state_array=ghost(state_array)

####

spl_psi = scipyinterpolate.splrep(initX[4:L-3], state_array[4:L-3,3],k=4)
psi_func(x) = scipyinterpolate.splev(x, spl_psi)

spl_derpsi = scipyinterpolate.splrep(initX[4:L-3], initderpsi[4:L-3],k=4)
derpsi_func(x) = scipyinterpolate.splev(x, spl_derpsi)

funcs = [psi_func derpsi_func]

#BETA
beta0=0
initbeta[4:L-3]=rk4wrapper(SFconstraint_beta,beta0,initX1,0,funcs)
state_array[:,2]=initbeta;
state_array=ghost(state_array);


#M
m0=0
initm[4:L-3]=rk4wrapper(SFconstraint_m,m0,initX1,0,funcs)
state_array[:,1]=initm;
state_array = ghost(state_array);

#CSV.write(dir*"/time_step0.csv", Tables.table(transpose(Matrix(state_array))), writeheader=false)
run=int(run)
CSV.write(dir*"/run$run/time_step0.csv", Tables.table(state_array), writeheader=false)

timestep=0
criticality=0.0
explode=0.0
critical_stop=0
evol_stats = [criticality A sigma r0 timestep explode run]
monitor_ratio = zeros(L)

run=int(run)
if run == 1
    CSV.write(dir*"/parameters.csv", Tables.table(evol_stats), writeheader=true, header=["criticality", "A", "sigma", "r0", "timestep", "explode", "run"])
end

finaltime=length(T)-1
stats=timeevolution(state_array,finaltime,dir,dt,run)


CSV.write("/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/DATA/bisectionsearch/parameters.csv", Tables.table(stats), writeheader=true,header=["criticality", "A", "sigma", "r0", "timestep", "explode", "run"],append=true);
