# Definition of gaussian initial data functions

function init_gaussian(x,r0,sigma,A)
    n=length(x);
    if n==1
        z= A * (x/(1-x))^2 * exp(-((x/(1-x)-r0)/sigma)^2)
    else
        z=zeros(n);
        for i in 1:n
            if i<n-4
                z[i] = A * (x[i]/(1-x[i]))^2 * exp(-((x[i]/(1-x[i])-r0)/sigma)^2)
            else
                z[n]=0 #avoid NaN for x=1
            end
        end
    end
    return z
end

function init_gaussian_der(x,r0,sigma,A)
    n=length(x);
    if n==1
        z= 2*A * x/(1-x)^3 * exp(-((x/(1-x)-r0)/sigma)^2) * (1 - x/(1-x) * (x/(1-x) - r0) / sigma^2)
    else
        z=zeros(n);
        for i in 1:n
            if i<n-4 #avoid NaN for x=1, otherwise, it's 0
                z[i] = 2*A * x[i]/(1-x[i])^3 * exp(-((x[i]/(1-x[i])-r0)/sigma)^2) * (1 - x[i]/(1-x[i]) * (x[i]/(1-x[i]) - r0) / sigma^2)
            end
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

# Updating Grid
function update_grid(data,dx,T,k,spl_funcs)
    
    #evolve grid
    data=rungekutta4molstep(Grid_RHS,data,T,k,0,spl_funcs) #evolve X
    data=ghost(data)

    X = data[:,5]

    return X
"""
    new_grid=[-3.0*dx, -2.0*dx, -dx]

    for i in X
        if i>=0 && i<=1
            new_grid = vcat(new_grid,i) #append
        end
    end
    new_grid = vcat(new_grid,[1.0+dx, 1.0+2.0*dx, 1.0+3.0*dx]) #append
    
    return new_grid"""
end
    
#Building initial data with a Runge-Kutta integrator for the constraint

function rungekutta4(f,y0,T)
    
    n = length(T)
    y = zeros(n)
    y[1] = y0;
    k_array=Array{Float64,2}(undef, 0, 5)
    for i in 1:n-1
        h = T[2] .- T[1]
        k1 = f(y[i], T[i])
        k2 = f(y[i] .+ k1 * h/2, T[i] .+ h/2)
        k3 = f(y[i] .+ k2 * h/2, T[i] .+ h/2)
        k4 = f(y[i] .+ k3 * h, T[i] .+ h)
        y[i+1] = y[i] .+ (h/6) * (k1 .+ 2*k2 .+ 2*k3 .+ k4)

        aux = [i k1 k2 k3 k4]
        k_array = vcat(k_array,aux)
    end
    return y,k_array
end



# Runge Kutta integrator used for the method of lines

function rungekutta4molstep(f,y00,T,w::Int64,ex,spl_func)
    X=state_array[:,5]
    y = y00;
        h = dt; # only for equally spaced grid in time and space, otherwise (T[w+1] - T[w])
        k1 = f(y[:,:], T[w],spl_func,X)
        k1=ghost(k1)
        k2 = f(y[:,:] + k1 * h/2, T[w] + h/2,spl_func,X)
        k2=ghost(k2)
        k3 = f(y[:,:] + k2 * h/2, T[w] + h/2,spl_func,X)
        k3=ghost(k3)
        k4 = f(y[:,:] + k3 * h, T[w] + h,spl_func,X)
        k4=ghost(k4)
        y[:,:] = y[:,:] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    return ghost(y[:,:])
end



"function rk4wrapper(f,y0,x,u) # u depicts T array or state_array data
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
end"

function rk4wrapper(f,y0,x,u,spl_func,psiarray) # u depicts T array, or M!!
    n = length(x)
    y = zeros(n)
    y[1] = y0;
    for i in 1:n-1
        h = x[2] .- x[1]
        k1 = f(y[i], x[i],u,spl_func,psiarray,2*i-1)
        k2 = f(y[i] .+ k1 * h/2, x[i] .+ h/2,u,spl_func,psiarray,2*i)
        k3 = f(y[i] .+ k2 * h/2, x[i] .+ h/2,u,spl_func,psiarray,2*i)
        k4 = f(y[i] .+ k3 * h, x[i] .+ h,u,spl_func,psiarray,2*i+1)
        y[i+1] = y[i] .+ (h/6) * (k1 .+ 2*k2 .+ 2*k3 .+ k4)
    end
    return y
end

#ghosts

function ghost(y)
    L=length(y[:,5])
    
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
function dissipation6(y,i,eps)
        delta6=(y[i+3,:]-6*y[i+2,:]+15*y[i+1,:]-20*y[i,:]+15*y[i-1,:]-6*y[i-2,:]+y[i-3,:]);
    return (-1)^3*eps*1/(dx)*delta6
end


#4th order  dissipation, added to 2nd order original scheme
function dissipation4(y,i,eps)
        delta4=(y[i+2,:]-4*y[i+1,:]+6*y[i,:]-4*y[i-1,:]+y[i-2,:]);
    return (-1)^2*eps*1/(dx)*delta4
end


#2nd order  dissipation, added to 1st order scheme ##NEW##
function dissipation2(y,i)
        delta2=(y[i+1,:]-2*y[i,:]+y[i-1,:]);
    return (-1)^1*epsilon*1/(dx)*delta2
end

# Discretization of derivatives
#Der(y,i,k,X)=(-y[i+2,k]+8*y[i+1,k]-8*y[i-1,k]+y[i-2,k])/(12*(X[i+1]-X[i])); #4th order
#DDer(y,i,k,X)=(-y[i+2,k]+16*y[i+1,k]-30*y[i,k]+16*y[i-1,k]-y[i-2,k])/(12*(X[i+1]-X[i])^2); #4th order

function Der(y,i,k,X)

    if i>=3 && i<=L-3
        z = (-y[i+2,k]+8*y[i+1,k]-8*y[i-1,k]+y[i-2,k])/(12*(X[i+1]-X[i]));
    else
        dx = X[2]-X[1]
        interp_data = []
        
        for j in 4:L-3
            aux = (-y[j+2,k]+8*y[j+1,k]-8*y[j-1,k]+y[j-2,k])/(12*(dx))
            interp_data = vcat(interp_data, aux)
        end
        spl = scipyinterpolate.splrep(X[4:L-3], interp_data,k=4)
        Der_func(x) = scipyinterpolate.splev(x, spl)
        z = Der_func(X[i])[1]
        #println("Der!, i=",i)
    end
    return z
end

function DDer(y,i,k,X)

    if i>=3 && i<=L-3
        z = (-y[i+2,k]+16*y[i+1,k]-30*y[i,k]+16*y[i-1,k]-y[i-2,k])/(12*(X[i+1]-X[i])^2);
    else
        dx = X[2]-X[1]
        interp_data = []
        
        for i in 4:L-3
            aux = (-y[i+2,k]+16*y[i+1,k]-30*y[i,k]+16*y[i-1,k]-y[i-2,k])/(12*(dx)^2)
            interp_data = vcat(interp_data, aux)
        end
        spl = scipyinterpolate.splrep(X[4:L-3], interp_data,k=4)
        DDer_func(xx) = scipyinterpolate.splev(xx, spl)
        z = DDer_func(X[i])[1]
        #println("DDer!, i=",i)
    end
    return z
end

function Der_cont(interp_func,x,i)

    dxx = X[i+1]-X[i]

    f=interp_func

    return ((-f(x+2*dxx) + 8*f(x+dxx) - 8*f(x-dxx) + f(x-2*dxx)) / (12*dxx))

end

function DDer_cont(interp_func,x,i)

    dxx = X[i+1]-X[i]

    f=interp_func

    return ((-f(x+2*dxx) + 16*f(x+dxx) - 30*f(x) + 16*f(x-dxx) - f(x-2*dxx)) / (12*dxx^2))

end


int(x) = floor(Int, x)


# Test Model RHSs for the bulk equations (3.6.16)


function SFconstraint_beta(beta0,x1,time,interp_func,psi_array,k)
    
    #psi = interp_funcs[3]
    derpsi = interp_func

    if x1<10^(-15)
        z = 0.0
    elseif abs.(x1 .- 1.0)<10^(-15)
        z = 0.0
    else
        #z = 2.0 .* pi .* (1.0 .- x1) ./ x1 .^3.0 .* (psi(x1) .+ (x1 .- 1) .* x1 .* derpsi(x1)) .^2
        z = 2.0 .* pi .* (1.0 .- x1) ./ x1 .^3.0 .* (psi_array[k] .+ (x1 .- 1) .* x1 .* derpsi(x1)) .^2
    end

    return z
end


function SFconstraint_m(m0,x1,u,interp_func,psiarray,k)

    derpsi = interp_func


    if x1<10^(-15)
        z = 0
    elseif abs.(x1 .- 1.0)<10^(-15)
        z = 2*pi .* (psiarray[k]) .^2
    else
        z = 2*pi .* (x1 .+ 2 .* (x1 .- 1) .* m0) ./ x1 .^3 .* (psiarray[k] .+ (x1 .- 1) .* x1 .* derpsi(x1)) .^2
    end
    
    return z
end
"
function SFconstraint_m(m0,x1,time,interp_funcs,new_m,i,bound)

    
    if i<=bound
        m=m = interp_funcs[1] #previous slice data
    else
        m=new_m #current slice data
    end

    psi = interp_funcs[3]
    derpsi = interp_funcs[4]


    if x1<10^(-15)
        z = 0
    elseif abs.(x1 .- 1.0)<10^(-15)
        z = 2*pi .* (psi(x1)) .^2
    else
        z = 2*pi .* (x1 .+ 2 .* (x1 .- 1) .* m(x1)) ./ x1 .^3 .* (psi(x1) .+ (x1 .- 1) .* x1 .* derpsi(x1)) .^2
    end
    
    return z
end"

function bulkSF(y,i,X)
    
    
    dy=zeros(length(y[1,:]));

    dy[1]=0; #m
    dy[2]=0; #beta
    dy[3]=0; #psi

    dy[4]=-1.0/2.0*exp(2.0*y[i,2])*((2.0*exp(2.0*(X[i]-y[i,3]+X[i]*y[i,3])*y[i,2]/X[i])*(X[i]-1.0)^2*(X[i]*((X[i]-1.0)*Der(y,i,1,X)+X[i]*Der(y,i,2,X))+y[i,1]*(1.0+2.0*(X[i]-1.0)*X[i]*Der(y,i,2,X))))/X[i]^2 - ((-1.0+X[i])^3*(X[i]+2.0*(X[i]-1.0)*y[i,1])*y[i,4])/X[i]^2 - ((1.0-X[i])^3*(1.0-2.0*(X[i]-1.0)^2*Der(y,i,1,X))*y[i,4])/X[i] - (2.0*(X[i]-1.0)^4*(X[i]+2.0*(X[i]-1.0)*y[i,1])*Der(y,i,2,X)*y[i,4])/X[i] - (Der(y,i,4,X)) - (2.0*(X[i]-1.0)*y[i,1]*Der(y,i,4,X))/X[i]);# - dissipation6(y,i,0.01)[4];
    
    dy[5]=0;

    return dy
end


function boundarySF(y,i)
    dy=zeros(length(y[1,:]));
    dy[1]=1.0; #f
    return dy
end

# Defining the function in the RHS of the evolution equation system

function SF_RHS(y,t,interp_func,X)
    
    L=length(X)
    dy=zeros((L,length(y[1,:])));

 
    for i in 5:(L-3)
        dy[i,4]=-1.0/2.0*exp(2.0*y[i,2])* ((2.0*exp(2.0*(X[i]-y[i,3]+X[i]*y[i,3])*y[i,2]/X[i])*(X[i]-1.0)^2*(X[i]*((X[i]-1.0)*Der(y,i,1,X)+X[i]*Der(y,i,2,X))+y[i,1]*(1.0+2.0*(X[i]-1.0)*X[i]*Der(y,i,2,X))))/X[i]^2 - ((-1.0+X[i])^3*(X[i]+2.0*(X[i]-1.0)*y[i,1])*y[i,4])/X[i]^2 - ((1.0-X[i])^3*(1.0-2.0*(X[i]-1.0)^2*Der(y,i,1,X))*y[i,4])/X[i] - (2.0*(X[i]-1.0)^4*(X[i]+2.0*(X[i]-1.0)*y[i,1])*Der(y,i,2,X)*y[i,4])/X[i] - (DDer(y,i,3,X)) - (2.0*(X[i]-1.0)*y[i,1]*DDer(y,i,3,X))/X[i]);# - dissipation6(y,i,0.01)[4];
    
        #termo a term
        #-1.0/2.0*exp(2.0*y[i,2])* (((-1.0+X[i])^3*(X[i]+2.0*(X[i]-1.0)*y[i,1])*y[i,4])/X[i]^2 - ((1.0-X[i])^3*(1.0-2.0*(X[i]-1.0)^2*Der(y,i,1,X))*y[i,4])/X[i] - (2.0*(X[i]-1.0)^4*(X[i]+2.0*(X[i]-1.0)*y[i,1])*Der(y,i,2,X)*y[i,4])/X[i] - (2.0*(X[i]-1.0)*y[i,1]*Der(y,i,4,X))/X[i]);

    end
    #expressao toda
    #-1.0/2.0*exp(2.0*y[i,2])* ((2.0*exp(2.0*(X[i]-y[i,3]+X[i]*y[i,3])*y[i,2]/X[i])*(X[i]-1.0)^2*(X[i]*((X[i]-1.0)*Der(y,i,1,X)+X[i]*Der(y,i,2,X))+y[i,1]*(1.0+2.0*(X[i]-1.0)*X[i]*Der(y,i,2,X))))/X[i]^2 - ((-1.0+X[i])^3*(X[i]+2.0*(X[i]-1.0)*y[i,1])*y[i,4])/X[i]^2 - ((1.0-X[i])^3*(1.0-2.0*(X[i]-1.0)^2*Der(y,i,1,X))*y[i,4])/X[i] - (2.0*(X[i]-1.0)^4*(X[i]+2.0*(X[i]-1.0)*y[i,1])*Der(y,i,2,X)*y[i,4])/X[i] - (Der(y,i,4,X)) - (2.0*(X[i]-1.0)*y[i,1]*Der(y,i,4,X))/X[i]);# - dissipation6(y,i,0.01)[4];
    

"""
    for i in 4:(L-3)
        if X[i]<10^(-7) #left #i<6

            #dy[i,1] to dy[i,3] stay 0
            dy[i,4]=0#-dissipation4(y,i,0.1)[4];

        elseif X[i] < (1-10^(-15)) #bulk
            dy[i,4]=-1.0/2.0*exp(2.0*y[i,2])*((2.0*exp(2.0*(X[i]-y[i,3]+X[i]*y[i,3])*y[i,2]/X[i])*(X[i]-1.0)^2*(X[i]*((X[i]-1.0)*Der(y,i,1,X)+X[i]*Der(y,i,2,X))+y[i,1]*(1.0+2.0*(X[i]-1.0)*X[i]*Der(y,i,2,X))))/X[i]^2 - ((-1.0+X[i])^3*(X[i]+2.0*(X[i]-1.0)*y[i,1])*y[i,4])/X[i]^2 - ((1.0-X[i])^3*(1.0-2.0*(X[i]-1.0)^2*Der(y,i,1,X))*y[i,4])/X[i] - (2.0*(X[i]-1.0)^4*(X[i]+2.0*(X[i]-1.0)*y[i,1])*Der(y,i,2,X)*y[i,4])/X[i] - (Der(y,i,4,X)) - (2.0*(X[i]-1.0)*y[i,1]*Der(y,i,4,X))/X[i]) - dissipation4(y,i,0.1)[4];

        else #right
            dy[i,4]=-1.0/2.0*exp(2.0*y[i,2])*((2.0*exp(2.0*(X[i]-y[i,3]+X[i]*y[i,3])*y[i,2]/X[i])*(X[i]-1.0)^2*(X[i]*((X[i]-1.0)*Der(y,i,1,X)+X[i]*Der(y,i,2,X))+y[i,1]*(1.0+2.0*(X[i]-1.0)*X[i]*Der(y,i,2,X))))/X[i]^2 - ((-1.0+X[i])^3*(X[i]+2.0*(X[i]-1.0)*y[i,1])*y[i,4])/X[i]^2 - ((1.0-X[i])^3*(1.0-2.0*(X[i]-1.0)^2*Der(y,i,1,X))*y[i,4])/X[i] - (2.0*(X[i]-1.0)^4*(X[i]+2.0*(X[i]-1.0)*y[i,1])*Der(y,i,2,X)*y[i,4])/X[i] - (Der(y,i,4,X)) - (2.0*(X[i]-1.0)*y[i,1]*Der(y,i,4,X))/X[i]) - dissipation4(y,i,0.1)[4];
        end
    end
    """
    
    #dy[4,:]=boundarySF(y,4);
    #dy[L-3,:]=boundarySF(y,L-3);

    #inner boundary
    
    
    #outer boundary
"""    dy[L-3,4]=extrapolate_out(dy[L-7,4], dy[L-6,4], dy[L-5,4], dy[L-4,4])#1.0/2.0*exp(2.0*y[L-3,2])* Der(y,L-3,4)
    dy[L-2,4]=extrapolate_out(dy[L-6,4], dy[L-5,4], dy[L-4,4], dy[L-3,4])#1.0/2.0*exp(2.0*y[L-3,2])* Der(y,L-2,4)
    dy[L-1,4]=extrapolate_out(dy[L-5,4], dy[L-4,4], dy[L-3,4], dy[L-2,4])
    dy[L,4]=extrapolate_out(dy[L-4,4], dy[L-3,4], dy[L-2,4], dy[L-1,4])"""
    

    return dy
end

function GP_RHS(y,t,interp_func,X)

    L=length(X)
    dy=zeros((L,length(y[1,:])));

    # updating ghostpoints
    for i in 1:3
        dy[i,4]= dy[i,4]=-1.0/2.0*exp(2.0*y[i,2])* ((2.0*exp(2.0*(X[i]-y[i,3]+X[i]*y[i,3])*y[i,2]/X[i])*(X[i]-1.0)^2*(X[i]*((X[i]-1.0)*Der(y,i,1,X)+X[i]*Der(y,i,2,X))+y[i,1]*(1.0+2.0*(X[i]-1.0)*X[i]*Der(y,i,2,X))))/X[i]^2 - ((-1.0+X[i])^3*(X[i]+2.0*(X[i]-1.0)*y[i,1])*y[i,4])/X[i]^2 - ((1.0-X[i])^3*(1.0-2.0*(X[i]-1.0)^2*Der(y,i,1,X))*y[i,4])/X[i] - (2.0*(X[i]-1.0)^4*(X[i]+2.0*(X[i]-1.0)*y[i,1])*Der(y,i,2,X)*y[i,4])/X[i] - (DDer(y,i,3,X)) - (2.0*(X[i]-1.0)*y[i,1]*DDer(y,i,3,X))/X[i]);# - dissipation6(y,i,0.01)[4];
    end
    for i in L-2:L
        dy[i,4]= dy[i,4]=-1.0/2.0*exp(2.0*y[i,2])* ((2.0*exp(2.0*(X[i]-y[i,3]+X[i]*y[i,3])*y[i,2]/X[i])*(X[i]-1.0)^2*(X[i]*((X[i]-1.0)*Der(y,i,1,X)+X[i]*Der(y,i,2,X))+y[i,1]*(1.0+2.0*(X[i]-1.0)*X[i]*Der(y,i,2,X))))/X[i]^2 - ((-1.0+X[i])^3*(X[i]+2.0*(X[i]-1.0)*y[i,1])*y[i,4])/X[i]^2 - ((1.0-X[i])^3*(1.0-2.0*(X[i]-1.0)^2*Der(y,i,1,X))*y[i,4])/X[i] - (2.0*(X[i]-1.0)^4*(X[i]+2.0*(X[i]-1.0)*y[i,1])*Der(y,i,2,X)*y[i,4])/X[i] - (DDer(y,i,3,X)) - (2.0*(X[i]-1.0)*y[i,1]*DDer(y,i,3,X))/X[i]);# - dissipation6(y,i,0.01)[4];
    end

    return dy
end

function Grid_RHS(y,t,interp_funcs,X)
    
    m = interp_funcs[1]
    beta = interp_funcs[2]

    L=length(X)
    dy=zeros((L,length(y[1,:])));

    

    for i in 4:(L-3)
        if X[i]<10^(-7) #left

            #dy[i,1] to dy[i,3] stay 0
            dy[i,5]=-1.0/2.0*exp(2.0*y[i,2]);

        else #bulk
            dy[i,5]=-1.0/2.0*(1-X[i])^2*exp(2.0*y[i,2])*(1-2*y[i,1]*(1-X[i])/X[i]);#dissipation4
        end
    end
    

    return dy
end

function doublegrid(X)
    new_grid=[X[1]]

    L = length(X)

    for i in 1:(L-1)
        h = X[i+1]-X[i]
        
        new_grid = vcat(new_grid, X[i]+h/2)
        new_grid = vcat(new_grid, X[i+1])
    end

    return new_grid

end