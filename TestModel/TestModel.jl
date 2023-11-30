# Definition of gaussian initial data functions

function scalar_timeder(R)
    n=length(R);
    if n==1
        z= (B* (exp(-((R - R01)/s1)^2) + exp(-((R + R01)/s1)^2))/(sqrt(2*pi)*s1))
    else
        z=zeros(n);
        for i in 1:n
            z[i]= (B* (exp(-((R[i] - R01)/s1)^2) + exp(-((R[i] + R01)/s1)^2))/(sqrt(2*pi)*s1))
        end
    end
    return z
end

function scalar_spaceder(R)
    n=length(R);
    if n==1
        z=(P* (-((2 *exp(-((R - R02)^2/s2^2)) *(R - R02))/s2^2) - (2*exp(-((R + R02)^2/s2^2))* (R + R02))/s2^2))/(sqrt(2*pi)* s2); #+ (P3* (-((2 *exp(-((R - R03)^2/s2^2)) *(R - R03))/s2^2) - (2*exp(-((R + R03)^2/s2^2))* (R + R03))/s2^2))/(sqrt(2*pi)* s2);
    else
    z=zeros(n);
    for i in 1:n
        z[i]=(P* (-((2 *exp(-((R[i] - R02)^2/s2^2)) *(R[i] - R02))/s2^2) - (2*exp(-((R[i] + R02)^2/s2^2))* (R[i] + R02))/s2^2))/(sqrt(2*pi)* s2); #+ (P3* (-((2 *exp(-((R[i] - R03)^2/s2^2)) *(R[i] - R03))/s2^2) - (2*exp(-((R[i] + R03)^2/s2^2))* (R[i] + R03))/s2^2))/(sqrt(2*pi)* s2)
    end
    end
    return z
end

function scalar_field(R)
    n=length(R);
    if n==1
        z= (P* (exp(-((R - R02)/s2)^2) + exp(-((R + R02)/s2)^2))/(sqrt(2*pi)*s2))
    else
        z=zeros(n);
        for i in 1:n
            z[i]= (P* (exp(-((R[i] - R02)/s2)^2) + exp(-((R[i] + R02)/s2)^2))/(sqrt(2*pi)*s2))
        end
    end
    return z
end


# Extrapolation

function extrapolate_in(y2,y3)
    return y2 + (y2-y3)
end

function extrapolate_out(y1,y2)
    return y2 + (y2-y1)
end


# NEW Extrapolation

function extrapolate_out_new(y0,y1,y2,y3)
    return -y0 + 4*y1 - 6*y2 + 4*y3
end

function extrapolate_in_new(y0,y1,y2,y3)
    return -y3 + 4*y2 - 6*y1 + 4*y0
end

    
#Building initial data with a Runge-Kutta integrator for the constraint

function rungekutta4(f,y0,T)
    
    n = length(T)
    y = zeros(n)
    y[1] = y0;
    for i in 1:n-1
        h = T[2] .- T[1]
        k1 = f(y[i], T[i])
        k2 = f(y[i] .+ k1 * h/2, T[i] .+ h/2)
        k3 = f(y[i] .+ k2 * h/2, T[i] .+ h/2)
        k4 = f(y[i] .+ k3 * h, T[i] .+ h)
        y[i+1] = y[i] .+ (h/6) * (k1 .+ 2*k2 .+ 2*k3 .+ k4)
    end
    return y
end


# Runge Kutta integrator used for the method of lines

function rungekutta4molstep(f,y00,T,w::Int64,ex)
    #y = zeros(length(R),2);
    y = y00;
        h = dt; # only for equally spaced grid in time and space, otherwise (T[w+1] - T[w])
        k1 = f(y[:,:], T[w])
        """println(T[w])
        println(T[w]+ h/2)
        println(T[w]+ h)"""
        k1=ghost(k1)
        k2 = f(y[:,:] + k1 * h/2, T[w] + h/2)
        k2=ghost(k2)
        k3 = f(y[:,:] + k2 * h/2, T[w] + h/2)
        k3=ghost(k3)
        k4 = f(y[:,:] + k3 * h, T[w] + h)
        k4=ghost(k4)
        y[:,:] = y[:,:] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    return ghost(y[:,:])
end

function rk4wrapper(f,y0,x,u)
    n = length(x)
    y = zeros(n)
    y[1] = y0;
    for i in 1:n-1
        h = x[2] .- x[1]
        k1 = f(y[i], x[i],u)
        k2 = f(y[i] .+ k1 * h/2, x[i] .+ h/2,u)
        k3 = f(y[i] .+ k2 * h/2, x[i] .+ h/2,u)
        k4 = f(y[i] .+ k3 * h, x[i] .+ h,u)
        y[i+1] = y[i] .+ (h/6) * (k1 .+ 2*k2 .+ 2*k3 .+ k4)
    end
    return y
end



using Printf
function print_muninn(files, t, data, res, mode, grid)
    #mode is "a" for append or "w" for write
    j=1
    #LL=length(state_array[1,:])
    if bisection==false
        for fl in files #normal run
            
            open(dir*"/$fl.txt", mode) do file
            #open("./DATA/muninnDATA/res$res/$fl.txt", mode) do file
                @printf file "\"Time = %.10e\n" t
                for i in 1:length(data[:,1])
                    @printf file "% .10e % .10e\n" grid[i] data[i,j]
                end
                println(file) # insert empty line to indicate end of data set
                end
            j=j+1
        end
    else
        if loggrid==true
            auxdir= dir*"/bisectionsearch/muninnDATA/uneven"
        else
            auxdir= dir*"/bisectionsearch/muninnDATA/even"
        end
        
        for fl in files #bisection search
            
            open(auxdir*"/run$run/$fl.txt", mode) do file
            #open("./DATA/bisectionsearch/muninnDATA/run$run/$fl.txt", mode) do file
                @printf file "\"Time = %.10e\n" t
                for i in 1:length(data[:,1])
                    @printf file "% .10e % .10e\n" grid[i] data[i,j]
                end
                println(file) # insert empty line to indicate end of data set
                end
            j=j+1
        end
    end
end



#ghosts

function ghost(y)
    L=length(y[:,1])
    
    #inner boundary extrapolation
    y[3,:]=extrapolate_in_new(y[4,:], y[5,:], y[6,:], y[7,:]) #new
    y[2,:]=extrapolate_in_new(y[3,:], y[4,:], y[5,:], y[6,:])
    y[1,:]=extrapolate_in_new(y[2,:], y[3,:], y[4,:], y[5,:])

    #outer boundary extrapolation
    y[L-2,:]=extrapolate_out_new(y[L-6,:], y[L-5,:], y[L-4,:], y[L-3,:]) #new
    y[L-1,:]=extrapolate_out_new(y[L-5,:], y[L-4,:], y[L-3,:], y[L-2,:])
    y[L,:]=extrapolate_out_new(y[L-4,:], y[L-3,:], y[L-2,:], y[L-1,:])
   

    return y
end


#6th order dissipation, added to 4th order original scheme #new
function dissipation6(y,i)
        delta6=(y[i+3,:]-6*y[i+2,:]+15*y[i+1,:]-20*y[i,:]+15*y[i-1,:]-6*y[i-2,:]+y[i-3,:]);
    return (-1)^3*epsilon*1/(dx)*delta6
end


#4th order  dissipation, added to 2nd order original scheme
function dissipation4(y,i)
        delta4=(y[i+2,:]-4*y[i+1,:]+6*y[i,:]-4*y[i-1,:]+y[i-2,:]);
    return (-1)^2*epsilon*1/(dx)*delta4
end


# Discretization of derivatives
#Der(y,i,k)=(y[i+1,k]-y[i-1,k])/(R[i+1]-R[i-1]); #2nd order
Der(y,i,k)=(-y[i+2,k]+8*y[i+1,k]-8*y[i-1,k]+y[i-2,k])/(12*(R[i+1]-R[i])); #4th order
#DDer(y,i,k)=(-y[i+2,k]+8*y[i+1,k]-8*y[i-1,k]-y[i-2,k])/(12*(R[i+1]-R[i])); #2nd order
DDer(y,i,k)=(-y[i+2,k]+16*y[i+1,k]-30*y[i,k]+16*y[i-1,k]-y[i-2,k])/(12*(R[i+1]-R[i])^2); #4th order


# Test Model RHSs for the bulk equations.
TMconstraint_f(f0,R1,time)=sin(R1 .+ time); ### function A(ut,xt)

function bulkTM(y,i)
    dy=zeros(length(y[1,:]));

    dy[1]=0;
    dy[2]=1.0/2.0*(0.0*y[i,1]+Der(y,i,2)) #g

    return dy
end

function bulkTMdiss(y,i)
    dy=zeros(length(y[1,:]));

    dy[1]=0;
    dy[2]=1.0/2.0*(0.0*y[i,1]+Der(y,i,2))-dissipation6(y,i)[2]; #g

    return dy
end

function boundaryTM(y,i)
    dy=zeros(length(y[1,:]));
    dy[1]=1.0; #f
    return dy
end

# Defining the function in the RHS of the evolution equation system

function TMRHS(y,t)
    L=length(R)
    dy=zeros((L,length(y[1,:])));


    y[4:L-3,1]=rk4wrapper(TMconstraint_f,f0,R1,t);
    
    
    y=ghost(y)

    global state_array[:,1] = y[:,1]
        
        for i in 4:(L-3)
                dy[i,:]=bulkTMdiss(y,i);
        end


    #dy[4,:]=boundaryTM(y,4);
    #dy[L-3,:]=boundaryTM(y,L-3);


    #outer boundary
    #dy[L-3,2]=2.0/10.0*pi*cos(2.0*pi/10.0*(R[L-2]*2.0+t));#y[L-2,1];
    dy[L-3,2]=2.0/10.0*pi*cos(2.0*pi/10.0*(R[L-3]*2.0+t))#+1/2*y[L-2,1]; # +1/2f, ie +1/2*y[l-2,1]
    dy[L-2,2]=extrapolate_out_new(dy[L-6,2], dy[L-5,2], dy[L-4,2], dy[L-3,2]) #new because I was having problems with g extrapolation at right border
    dy[L-1,2]=extrapolate_out_new(dy[L-5,2], dy[L-4,2], dy[L-3,2], dy[L-2,2])
    dy[L,2]=extrapolate_out_new(dy[L-4,2], dy[L-3,2], dy[L-2,2], dy[L-1,2])
    return dy
end

