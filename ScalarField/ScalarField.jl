# Definition of gaussian initial data functions

function scalar_spaceder(R)
    n=length(R);
    if n==1
        z=(Amp* (-((2 *exp(-((R - c)^2/sigma^2)) *(R - c))/sigma^2) - (2*exp(-((R + c)^2/sigma^2))* (R + c))/sigma^2))/(sqrt(2*pi)* sigma);
    else
    z=zeros(n);
    for i in 1:n
        z[i]=(Amp* (-((2 *exp(-((R[i] - c)^2/sigma^2)) *(R[i] - c))/sigma^2) - (2*exp(-((R[i] + c)^2/sigma^2))* (R[i] + c))/sigma^2))/(sqrt(2*pi)* sigma);
    end
    end
    return z
end

function scalar_field(R)
    n=length(R);
    if n==1
        z= (Amp* (exp(-((R - c)/sigma)^2) + exp(-((R + c)/sigma)^2))/(sqrt(2*pi)*sigma))
    else
        z=zeros(n);
        for i in 1:n
            z[i]= (Amp* (exp(-((R[i] - c)/sigma)^2) + exp(-((R[i] + c)/sigma)^2))/(sqrt(2*pi)*sigma))
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
        k1 = f(y[:,:], T[w],y00)
        """println(T[w])
        println(T[w]+ h/2)
        println(T[w]+ h)"""
        k1=ghost(k1)
        k2 = f(y[:,:] + k1 * h/2, T[w] + h/2,y00)
        k2=ghost(k2)
        k3 = f(y[:,:] + k2 * h/2, T[w] + h/2,y00)
        k3=ghost(k3)
        k4 = f(y[:,:] + k3 * h, T[w] + h,y00)
        k4=ghost(k4)
        y[:,:] = y[:,:] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    return ghost(y[:,:])
end

"""function rk4wrapper(f,y0,x,u, statearray_data)
    n = length(x)
    y = zeros(n)
    y[1] = y0;
    for i in 1:n-1
        h = x[2] .- x[1]
        k1 = f(y[i], x[i],u, statearray_data)
        k2 = f(y[i] .+ k1 * h/2, x[i] .+ h/2,u, statearray_data)
        k3 = f(y[i] .+ k2 * h/2, x[i] .+ h/2,u, statearray_data)
        k4 = f(y[i] .+ k3 * h, x[i] .+ h,u, statearray_data)
        y[i+1] = y[i] .+ (h/6) * (k1 .+ 2*k2 .+ 2*k3 .+ k4)
    end
    return y
end"""

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


# Test Model RHSs for the bulk equations (3.6.16)

function SFconstraint_beta(beta0,R1,time) # y is statearray_data
    if R1<10^(-7)
        z = 0
        print("R=0!!!")
    else
        #z = 2 .* pi .* (1 .- R1) ./ R1 .^3 .* (state_array[int.(R1 ./ dx .+ 1),3] .+ (R1 .- 1) .* R1 .* state_array[int.(R1 ./ dx .+ 1),4]) .^2
        z = 2 .* pi .* (1 ./ R1.^3 .- 1 ./ R1.^2) .* (state_array[int.(R1 ./ dx .+ 1),3] .+ (R1 .- 1) .* R1 .* state_array[int.(R1 ./ dx .+ 1),4]) .^2
        #z = 2 .* pi .* (1 .- R1) ./ R1 .^3 .* (y[int.(R1 ./ dx .+ 4),3] .+ (R1 .- 1) .* R1 .* y[int.(R1 ./ dx .+ 4),4]) .^2
    end


    #print(y[:,1])
    return z
end


function SFconstraint_m(m0,R1,time)
    if R1<10^(-7)
        z = 0
        print("R=0!!!")
    else
        z = 2*pi .* (R1 .+ 2 .* (R1 .- 1) .* state_array[int.(R1 ./ dx .+ 1),1]) ./ R1 .^3 .* (state_array[int.(R1 ./ dx .+ 1),3] .+ (R1 .- 1) .* R1 .* state_array[int.(R1 ./ dx .+ 1),4]) .^2
        #z = 2*pi .* (R1 .+ 2 .* (R1 .- 1) .* y[int.(R1 ./ dx .+ 4),1]) ./ R1 .^3 .* (y[int.(R1 ./ dx .+ 4),3] .+ (R1 .- 1) .* R1 .* y[int.(R1 ./ dx .+ 4),4]) .^2
    end

    return z
end

"function globalfunc()
    return  R1 .* state_array[4:L-3,4]
end"

function bulkSF(y,i)
    dy=zeros(length(y[1,:]));

    dy[1]=0; #m
    dy[2]=0; #beta
    dy[3]=0; #psi
    
    if i<5 #for i<3 I get a NaN
        dy[4]=5  # WRITE BOXOPERATOR APPROX HERE

    elseif R[i]>0.95
        #dy[4]=1.0/2.0*exp(2.0*y[i,2])* DDer(y,i,3); #psi,x
        dy[4]=5


    else
        #dy[4]=1.0/2.0*exp(2.0*y[i,2])* DDer(y,i,3); #psi,x
        dy[4]=5


        """println("i = ", i, " dy[4] = ", dy[4])
        println(" DDer(y,i,3) = ", DDer(y,i,3))
        println("num = ", (-y[i+2,3]+16*y[i+1,3]-30*y[i,3]+16*y[i-1,3]-y[i-2,3]))
        println(-y[i+2,3],16*y[i+1,3],30*y[i,3],16*y[i-1,3],y[i-2,3])"""
        

        #dy[4]=-1.0/2.0*exp(2.0*y[i,2])*((2.0*exp(2.0*(R[i]-y[i,3]+R[i]*y[i,3])*y[i,2]/R[i])*(R[i]-1.0)^2*(R[i]*((R[i]-1.0)*Der(y,i,1)+R[i]*Der(y,i,2))+y[i,1]*(1.0+2.0*(R[i]-1.0)*R[i]*Der(y,i,2))))/R[i]^2 - ((-1.0+R[i])^3*(R[i]+2.0*(R[i]-1.0)*y[i,1])*Der(y,i,3))/R[i]^2 - ((1.0-R[i])^3*(1.0-2.0*(R[i]-1.0)^2*Der(y,i,1))*Der(y,i,3))/R[i] - (2.0*(R[i]-1.0)^4*(R[i]+2.0*(R[i]-1.0)*y[i,1])*Der(y,i,2)*Der(y,i,3))/R[i] - (DDer(y,i,3)) - (2.0*(R[i]-1.0)*y[i,1]*DDer(y,i,3))/R[i]); #psi,x
        #dy[4]=1.0/2.0*exp(2.0*y[i,2])* DDer(y,i,3); #psi,x
    end
    
    return dy
end

function bulkSFdiss(y,i)
    dy=zeros(length(y[1,:]));

    dy[1]=0; #m
    dy[2]=0; #beta
    dy[3]=0; #psi

    #dy[4]=-1.0/2.0*exp(2.0*beta)*((2.0*exp(2*(x-psi+x*psi)*beta/x)*(x-1)^2*(x*((x-1)*Dm+x*Dbeta)+m*(1+2*(x-1)*x*Dbeta)))/x^2 - ((-1.0+x)^3*(x+2.0*(x-1.0)*m)*Dpsi)/x^2 - ((1-x)^3*(1-2*(x-1)^2*Dm)*Dpsi)/x - (2*(x-1)^4*(x+2(x-1)*m)*Dbeta*Dpsi)/x - (DDpsi)- (2*(x-1)*m*DDpsi)/x)
    # x = R[i]
    # m = y[i,1]
    # beta = y[i,2]
    # Dm = Der(y,i,1)
    # Dbeta = Der(y,i,2)
    # Dpsi = Der(y,i,3)
    # DDpsi = DDer(y,i,3)
    
    if i>3
        dy[4]=-1.0/2.0*exp(2.0*y[i,2])*((2.0*exp(2.0*(R[i]-y[i,3]+R[i]*y[i,3])*y[i,2]/R[i])*(R[i]-1.0)^2*(R[i]*((R[i]-1.0)*Der(y,i,1)+R[i]*Der(y,i,2))+y[i,1]*(1.0+2.0*(R[i]-1.0)*R[i]*Der(y,i,2))))/R[i]^2 -
        ((-1.0+R[i])^3*(R[i]+2.0*(R[i]-1.0)*y[i,1])*Der(y,i,3))/R[i]^2 - ((1.0-R[i])^3*(1.0-2.0*(R[i]-1.0)^2*Der(y,i,1))*Der(y,i,3))/R[i] - (2.0*(R[i]-1.0)^4*(R[i]+2.0*(R[i]-1.0)*y[i,1])*Der(y,i,2)*Der(y,i,3))/R[i]
        - (DDer(y,i,3)) - (2.0*(R[i]-1.0)*y[i,1]*DDer(y,i,3))/R[i]) - dissipation6(y,i)[4]; #psi,x
    else
        dy[4]=1 # WRITE BOXOPERATOR APPROX HERE
    end

    return dy
end

function boundarySF(y,i)
    dy=zeros(length(y[1,:]));
    dy[1]=1.0; #f
    return dy
end

# Defining the function in the RHS of the evolution equation system

function SF_RHS(y,t, statearray_data)
    L=length(R)
    dy=zeros((L,length(y[1,:])));


    #y[4:L-3,1]=rk4wrapper(SFconstraint_f,f0,R1,t,statearray_data);    
    #y=ghost(y)

    #global state_array[:,1] = y[:,1]
        
    for i in 4:(L-3)
            dy[i,:]=bulkSF(y,i);
            #dy[i,4]=bulkSFdiss(y,i);
    end


    #dy[4,:]=boundarySF(y,4);
    #dy[L-3,:]=boundarySF(y,L-3);

    #inner boundary
    dy[1,4]=0
    dy[2,4]=0#try dy[2,:]=0
    dy[3,4]=0
    dy[4,4]=0
    
    #outer boundary
    dy[L-3,4]=0#?
    dy[L-2,4]=0
    dy[L-1,4]=0
    dy[L,4]=0
    
    #dy[L-3,2]=2.0/10.0*pi*cos(2.0*pi/10.0*(R[L-2]*2.0+t));#y[L-2,1];
    "dy[L-3,2]=2.0/10.0*pi*cos(2.0*pi/10.0*(R[L-3]*2.0+t))#+1/2*y[L-2,1]; # +1/2f, ie +1/2*y[l-2,1]
    dy[L-2,2]=extrapolate_out(dy[L-6,2], dy[L-5,2], dy[L-4,2], dy[L-3,2])
    dy[L-1,2]=extrapolate_out(dy[L-5,2], dy[L-4,2], dy[L-3,2], dy[L-2,2])
    dy[L,2]=extrapolate_out(dy[L-4,2], dy[L-3,2], dy[L-2,2], dy[L-1,2])"
    return dy
end

