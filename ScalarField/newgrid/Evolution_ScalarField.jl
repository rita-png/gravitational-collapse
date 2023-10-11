using CSV, Tables
using Dierckx


A = ARGS[1]
run = trunc(Int, ARGS[2])
N = ARGS[3]


m = 1
res=m


function compactify(r)
    return r/(1+r)
end

function uncompactify(x)
    return x/(1-x)
end


if compactified==true
    Xf=1.0
else
    Xf=10.0
end

dx=Xf/N#Float128(Xf/N);
if loggrid==false
    dt=0.5*round(dx,digits=10)#0.5*dx#round(dx,digits=10);#dx
else
    dt=0.1*round(dx,digits=10)
end
Nt=N
Tf=Nt*dt; #final time
#print(Tf)

using Printf

if loggrid==true
    ori=(tan(-mgrid/Agrid)/fgrid+kgrid)#0.0#Float128(0.0)#0.0;
    Xf=(tan((1-mgrid)/Agrid)/fgrid+kgrid)
else
    ori=0.0
    Xf=1.0
end

dx=(Xf-ori)/N

initX1 = nothing

initX1=range(ori, stop=Xf, step=dx);
#initX1=create_range(ori,Xf,dx,N)
#initX = range(round(ori-3.0*dx,digits=10), stop=Xf+3.0*dx, step=dx)
#initX=create_range(ori-3.0*dx,Xf+3.0*dx,dx,N+6)

L=length(initX1)+6;#length(initX)

initX=[ori-3*dx; ori-2*dx; ori-dx; collect(initX1); Xf+dx; Xf+2*dx; Xf+3*dx];

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
if loggrid==true
    derpsi_func = Spline1D(initX[4:L-3], state_array[4:L-3,4],  k=4);
else
    derpsi_func = Spline1D(initX[4:L-3], state_array[4:L-3,4],  k=4);
end;


y0=[0.0, 0.0, 0.0]

if twod==false
    state_array[4:L-3,1:3] = n_rk4wrapper(RHS,y0,initX[4:L-3],0,derpsi_func,state_array[:,:]);
else
    state_array[4:L-3,1:3] = twod_n_rk4wrapper(RHS,y0,initX[4:L-3],0,derpsi_func,state_array[:,:]);
end


state_array = ghost(state_array);

#threshold for apparent black hole formation
global monitor_ratio=zeros(L)
if compactified==false
    global monitor_ratio[5:L-4] = 2 .* state_array[5:L-4,1] ./ initX[5:L-4]
else
    if loggrid==false
        global monitor_ratio[5:L-4] = 2 .* state_array[5:L-4,1] ./ initX[5:L-4] .* (1 .- initX[5:L-4])
    else
        global monitor_ratio[5:L-4] = 2 .* state_array[5:L-4,1] ./ inverse(initX[5:L-4]) .* (1 .- inverse(initX[5:L-4]))
    end
end

global files=["m", "beta", "psi", "derpsi", "monitorratio"]


if zeroformat==true
    zero_print_muninn(files, 0, state_array[:,1:5],res,"w")
else
    print_muninn(files, 0, state_array[:,1:5],res,"w")
end


time=0.0
criticality=0.0
explode=0.0
maxever_mr=0.0
evol_stats = [criticality A sigma r0 time explode run maxever_mr]
global monitor_ratio = zeros(L);

run=int(run)
if run == 1 && bisection==true
    if loggrid==true
        CSV.write(dir*"/bisectionsearch/muninnDATA/uneven/parameters.csv", Tables.table(evol_stats))#, writeheader=true, header=["criticality", "A", "sigma", "r0", "time", "explode", "run"])
    else
        CSV.write(dir*"/bisectionsearch/muninnDATA/even/parameters.csv", Tables.table(evol_stats))#, writeheader=true, header=["criticality", "A", "sigma", "r0", "time", "explode", "run"])
    end
end

ginit=speed(initX,state_array[:,1],state_array[:,2])

using Base.Threads
Threads.nthreads()

#global dt=5e-5

if bisection==false
    if m==1
        #global dt=2e-5/5
        global dt=0.000006 #this is not the file  that is read for convergence tests
    elseif m==2
        #global dt=2e-5/5/2
        global dt=0.000006/2
    else
        #global dt=2e-5/5/2/2
        global dt=0.000006/2/2
    end
    finaltime=5.0
else
    finaltime=2.2
end

stats, T_interp = timeevolution(state_array,finaltime,run);

if bisection==true
    if loggrid==true
        CSV.write(dir*"/bisectionsearch/muninnDATA/uneven/parameters.csv", Tables.table(stats),append=true);#, writeheader=true,header=["criticality", "A", "sigma", "r0", "time", "explode", "run"],
        CSV.write(dir*"/bisectionsearch/muninnDATA/uneven/timearray.csv", Tables.table(T_interp))#, writeheader=false);
    else
        CSV.write(dir*"/bisectionsearch/muninnDATA/even/parameters.csv", Tables.table(stats),append=true)#, writeheader=true,header=["criticality", "A", "sigma", "r0", "time", "explode", "run"]);
        CSV.write(dir*"/bisectionsearch/muninnDATA/even/timearray.csv", Tables.table(T_interp))#, writeheader=false);
    end

    
end