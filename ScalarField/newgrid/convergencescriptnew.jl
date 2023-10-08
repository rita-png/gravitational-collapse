A = ARGS[1]
res = trunc(Int, ARGS[2])
#A=0.1
#m=1
#res=m
global compactified=true
global loggrid=true
global zeroformat = false
global twod = false
global bisection=false
m=res
run = 1


function compactify(r)
    return r/(1+r)
end

function uncompactify(x)
    return x/(1-x)
end

if compactified==true
    Xf=1.0
else
    Xf=10.0#Float128(1.0);
end

using Printf


res=m;
N=2.0^m*100.0/2.0
ori=0.0
Xf=1.0
dx=(Xf-ori)/N

println("running for resolution ", res, " N1 = ", N, ", A = ", A)

global dir = "/home/ritapsantos/data/ritapsantos/chebyconvergence"

using Printf
include("./ScalarField.jl");
#include("/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/myspline.jl");

ori=0.0#Float128(0.0)#0.0;
initX1 = nothing
N=int(N)
initX1=range(ori, stop=Xf, step=dx);
#initX1=create_range(ori,Xf,dx,N)
initX = range(round(ori-3.0*dx,digits=10), stop=Xf+3.0*dx, step=dx)
#initX=create_range(ori-3.0*dx,Xf+3.0*dx,dx,N+6)

L=length(initX);
println("step size is  ", dx)

using Dierckx
if loggrid==true
    global originalX=initX
    xtilde=gridfunc(initX1)
    initX1=xtilde
    initX=collect(initX)
    initX[4:L-3]=xtilde
    #global dergrid_func = der_grid(initX)
    global jacobian_func = Spline1D(originalX[4:L-3], analytic_jacobian(originalX[4:L-3]),  k=4);
end;
#initX=[ori-3*dx; ori-2*dx; ori-dx; collect(initX1); Xf+dx; Xf+2*dx; Xf+3*dx];




####

initm=zeros(L)
initbeta=zeros(L)
initpsi=zeros(L)
initderpsi=zeros(L)

state_array=[initm initbeta initpsi initderpsi initX];

#PSI
r0=0.7#Float128(0.7)#0.01#0.7#0.01#0.7#0.7#0.7#0.01#0.7#0.3
sigma=0.3#Float128(0.3)


#PSI,R FROM PSI
initderpsi[4:L-3] = init_gaussian_der(initX1,r0,sigma,A)


state_array[:,4] = initderpsi
state_array=ghost(state_array)

####
"""if loggrid==true
    derpsi_func = Spline1D(inverse(initX[4:L-3]), state_array[4:L-3,4],  k=4);
else
    derpsi_func = Spline1D(initX[4:L-3], state_array[4:L-3,4],  k=4);
end;"""

derpsi_func = Spline1D(initX[4:L-3], state_array[4:L-3,4],  k=4);

y0=[0.0 0.0 0.0]

if twod==false
    state_array[4:L-3,1:3] = n_rk4wrapper(RHS,y0,initX[4:L-3],0,derpsi_func,state_array[:,:]);
else
    state_array[4:L-3,1:3] = twod_n_rk4wrapper(RHS,y0,initX[4:L-3],0,derpsi_func,state_array[:,:]);
end

state_array = ghost(state_array);





using Tables

global files=["m", "beta", "psi", "derpsi"]

print_muninn(files, 0, state_array[:,1:5],res,"w")

time=0.0
criticality=0.0
explode=0.0
evol_stats = [criticality A sigma r0 time explode run]
global monitor_ratio = zeros(L);
#CSV.write(dir*"/parameters.csv", Tables.table(evol_stats), writeheader=true, header=["criticality", "A", "sigma", "r0", "time", "explode", "run"])

ginit=speed(initX,state_array[:,1],state_array[:,2])


using Base.Threads
Threads.nthreads()


if bisection==false
    if m==1
        #global dt=0.0004 #N=100
        global dt=5e-5 #N=100
    elseif m==2
        #global dt=0.0002 #N=200
        global dt=5e-5/2
    elseif m==3
        #global dt=0.0002/2 #N=200
        global dt=5e-5/2/2
    elseif m==4
        #global dt=0.0002/2/2 #N=200
        global dt=5e-5/2/2/2
    else
        #global dt=0.0002/2/2/2 #N=1600
        global dt=5e-5/2/2/2/2
    end
    finaltime=2.5
else
    finaltime=2.5
end


evol_stats, T_interp = timeevolution(state_array,finaltime,run);