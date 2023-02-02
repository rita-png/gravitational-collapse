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

# Interpolation

function interpolate(x,x1,x2,y1,y2)
    return y1+(y2-y1)*(x-x1)/(x2-x1)
end
#(1−x)^4=x^4−4x^3+6x^2−4x+1 ex. 4th order polynomial



# Extrapolation

function extrapolate_out(y0,y1,y2,y3)
    return -y0 + 4*y1 - 6*y2 + 4*y3
end

function extrapolate_in(y0,y1,y2,y3)
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


"function rungekutta4_data(fbar_data,y0,x)
    
    n = length(x)
    y = zeros(n)
    y[1] = y0;
    for i in 1:n-1
        h = x[2] .- x[1]
        k1 = f(y[i], x[i])
        k1 = fbar_data[i]

        k2 = f(y[i] .+ k1 * h/2, x[i] .+ h/2)
        k2 = interpolate(x[i],x[i+1],fbar_data[i],fbar_data[i+1])

        k3 = f(y[i] .+ k2 * h/2, x[i] .+ h/2)
        k4 = f(y[i] .+ k3 * h, x[i] .+ h)
        y[i+1] = y[i] .+ (h/6) * (k1 .+ 2*k2 .+ 2*k3 .+ k4)
    end
    return y
end"


#ghosts

function ghost(y)
    L=length(y[:,1])
    
    #inner boundary extrapolation
    y[3,:]=extrapolate_in(y[4,:], y[5,:], y[6,:], y[7,:])
    y[2,:]=extrapolate_in(y[3,:], y[4,:], y[5,:], y[6,:])
    y[1,:]=extrapolate_in(y[2,:], y[3,:], y[4,:], y[5,:])

    #outer boundary extrapolation
    y[L-2,:]=extrapolate_out(y[L-6,:], y[L-5,:], y[L-4,:], y[L-3,:])
    y[L-1,:]=extrapolate_out(y[L-5,:], y[L-4,:], y[L-3,:], y[L-2,:])
    y[L,:]=extrapolate_out(y[L-4,:], y[L-3,:], y[L-2,:], y[L-1,:])
   

    return y
end


#6th order dissipation, added to 4th order original scheme
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
Der(y,i,k)=(-y[i+2,k]+8*y[i+1,k]-8*y[i-1,k]+y[i-2,k])/(12*(R[i+1]-R[i])); #4th order
DDer(y,i,k)=(-y[i+2,k]+16*y[i+1,k]-30*y[i,k]+16*y[i-1,k]-y[i-2,k])/(12*(R[i+1]-R[i])^2); #4th order

int(x) = floor(Int, x)

# Initial data for psibar and psi
SFconstraint_psibar(psibar0,R1)=-sin.(R1)

# Test Model RHSs for the bulk equations.
# (3.6.16)
#SFconstraint_m(m0,R1,time)=2*pi .* R1 .^2 *(1 .- 2 .* (1 .- R1) ./ R1 .* state_array[int.(R1 ./ dx .+ 1),1]) * ((1 .- R1) ./ R1 .* state_array[int.(R1./dx.+1),4] .- state_array[int.(R1./dx.+1),3] ./ R1^2)^2

SFconstraint_m(beta0,R1,time)=2*pi ./ R1 .^2 .* (1-2*(1 .- R1) .* state_array[int.(R1 ./ dx .+ 1),1] ./ R1) .* (state_array[int.(R1 ./ dx .+ 1),3] .+ (R1 .- 1) R1 .* state_array[int.(R1 ./ dx .+ 1),4]) .^2

#SFconstraint_beta(beta0,R1,time)=2*pi .* R1 .* (1 .- R1) * ((1 .- R1) ./ R1 .* state_array[int.(R1./dx.+1),4] .- state_array[int.(R1./dx.+1),3] ./ R1^2)^2
SFconstraint_beta(m0,R1,time)=2*pi ./ R1 .^3 .*(R1 .+ 2 .* (R1 .- 1) .* state_array[int.(R1 ./ dx .+ 1),1] ) .* (state_array[int.(R1 ./ dx .+ 1),3] .+ (R1 .- 1) .* R1 .* state_array[int.(R1 ./ dx .+ 1),4]) .^2


"function globalfunc()
    return  R1 .* state_array[4:L-3,4]
end"

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
    dy[L-2,2]=extrapolate_out(dy[L-6,2], dy[L-5,2], dy[L-4,2], dy[L-3,2])
    dy[L-1,2]=extrapolate_out(dy[L-5,2], dy[L-4,2], dy[L-3,2], dy[L-2,2])
    dy[L,2]=extrapolate_out(dy[L-4,2], dy[L-3,2], dy[L-2,2], dy[L-1,2])
    return dy
end

