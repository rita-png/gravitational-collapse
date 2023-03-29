# Definition of gaussian initial data functions

function init_gaussian(x,r0,sigma,A)
    n=length(x);
    if n==1
        z= A * (x/(1-x))^3.0 * exp(-((x/(1-x)-r0)/sigma)^2)
    else
        z=zeros(n);
        for i in 1:n
            if i<n-4
                z[i] = A * (x[i]/(1-x[i]))^3.0 * exp(-((x[i]/(1-x[i])-r0)/sigma)^2.0)
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
        #z= 2*A * x/(1-x)^3 * exp(-((x/(1-x)-r0)/sigma)^2) * (1 - x/(1-x) * (x/(1-x) - r0) / sigma^2)
        z= A * exp(-((x/(1-x)-r0)/sigma)^2) * (3 * x^2 / (1-x) ^4 - (x/(1-x))^3 * (2*(x-r0*(-x+1)))/(sigma^2*(1-x)^3))
    else
        z=zeros(n);
        for i in 1:n
            if i<n-4 #avoid NaN for x=1, otherwise, it's 0
                #z[i] = 2*A * x[i]/(1-x[i])^3 * exp(-((x[i]/(1-x[i])-r0)/sigma)^2) * (1 - x[i]/(1-x[i]) * (x[i]/(1-x[i]) - r0) / sigma^2)
                z[i]=A * exp(-((x[i]/(1-x[i])-r0)/sigma)^2) * (3 * x[i]^2 / (1-x[i]) ^4 - (x[i]/(1-x[i]))^3 * (2*(x[i]-r0*(-x[i]+1)))/(sigma^2*(1-x[i])^3))
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
function update_grid(data,T,k)
    
    dx = data[2,5]-data[1,5]

    #evolve grid
    data = rungekutta4molstep(Grid_RHS,data,T,k,0,data[:,5]) #evolve X
    #data=ghost(data)
    
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

function rungekutta4molstep(f,y00,T,w::Int64,ex,X)
    y = y00;
    X=collect(X)
        h = T[w+1]-T[w]
        #print("\n\nh = ", h, " \n")
        #h = dt; # only for equally spaced grid in time and space, otherwise (T[w+1] - T[w])
        #println("y[:,4] ", y[:,4])
        #println("############")
        k1 = f(y[:,:], T[w],X)
        #k1=ghost(k1)
        k1=boundarySF(k1,X)
        #println("k1[:,4] ", k1[:,4])
        ####println("y[:,4] ", y[:,4])
        #println("############")
        k2 = f(y[:,:] .+ k1 .* h/2, T[w] + h/2,X)
        k2=boundarySF(k2,X)
        #println("k2[:,4] ", k2[:,4])
        ####println("y[:,4] .+ k1 .* h/2 ", y[:,4] .+ k1 .* h/2)
        #println("############")

        #println("k2[:,4].-k1[:,4] ", k2[:,4] .- k1[:,4])
        ####println("y[:,4] .+ k1 .* h/2 .- y[:,4]", y[:,4] .+ k1 .* h/2 .- y[:,4])
        #println("############")
        k3 = f(y[:,:] .+ k2 .* h/2, T[w] + h/2,X)
        k3=boundarySF(k3,X)

        #println("############")
        k4 = f(y[:,:] .+ k3 .* h, T[w] + h,X)
        k4=boundarySF(k4,X)

        """print("\n\ny[50,4] .+ k3 .* h\n\n",y[50,4] .+ k3[50,4] .* h)
        print("\n\ny[50,4]\n\n",y[50,4])
        print("\n\nk3\n\n",k3[50,4])
        print("\n\nh\n\n",h)
"""
        #println("k4[:,4] .- k1[:,4] ", k4[:,4] .- k1[:,4])
        #println("############")

        #println(y[:,4]-.y[:,4] .- k1 .* h/2)
        y[:,:] = y[:,:] .+ (h/6) .* (k1 .+ 2 * k2 .+ 2 * k3 .+ k4)
        
    return y[:,:]#ghost(y[:,:])
end



function rk4wrapper(f,y0,x,u,spl_funcs) # u depicts T array, or M!!
    n = length(x)
    y = zeros(n)
    y[1] = y0;
    for i in 1:n-1
        h = x[2] .- x[1]
        k1 = f(y[i], x[i],u,spl_funcs)
        k2 = f(y[i] .+ k1 * h/2, x[i] .+ h/2,u,spl_funcs)
        k3 = f(y[i] .+ k2 * h/2, x[i] .+ h/2,u,spl_funcs)
        k4 = f(y[i] .+ k3 * h, x[i] .+ h,u,spl_funcs)
        y[i+1] = y[i] .+ (h/6) * (k1 .+ 2*k2 .+ 2*k3 .+ k4)
    end
    return y
end

#ghosts

function ghost(y)
    L=length(y[:,1])
    
    #inner boundary extrapolation
    y[3,1:4]=extrapolate_in(y[4,1:4], y[5,1:4], y[6,1:4], y[7,1:4])
    y[2,1:4]=extrapolate_in(y[3,1:4], y[4,1:4], y[5,1:4], y[6,1:4])
    y[1,1:4]=extrapolate_in(y[2,1:4], y[3,1:4], y[4,1:4], y[5,1:4])

    #outer boundary extrapolation
    y[L-2,1:4]=extrapolate_out(y[L-6,1:4], y[L-5,1:4], y[L-4,1:4], y[L-3,1:4])
    y[L-1,1:4]=extrapolate_out(y[L-5,1:4], y[L-4,1:4], y[L-3,1:4], y[L-2,1:4])
    y[L,1:4]=extrapolate_out(y[L-4,1:4], y[L-3,1:4], y[L-2,1:4], y[L-1,1:4])
   

    return y
end


#6th order dissipation, added to 4th order original scheme
"""function dissipation6(y,i,eps)
        delta6=(y[i+3,:]-6*y[i+2,:]+15*y[i+1,:]-20*y[i,:]+15*y[i-1,:]-6*y[i-2,:]+y[i-3,:]);
    return (-1)^3*eps*1/(dx)*delta6
end"""

function dissipation6(y,i,eps)
    if i==4
        delta6= (19*y[i,:]-142*y[i+1,:]+464*y[i+2,:]-866*y[i+3,:]+1010*y[i+4,:]-754*y[i+5,:]+352*y[i+6,:]-94*y[i+7,:]+11*y[i+8,:])/2;
    elseif i==5
        delta6= (11*y[i-1,:]-80*y[i,:]+254*y[i+1,:]-460*y[i+2,:]+520*y[i+3,:]-376*y[i+4,:]+170*y[i+5,:]-44*y[i+6,:]+5*y[i+7,:])/2;
    elseif i==6
        delta6= (5*y[i-2,:]-34*y[i-1,:]+100*y[i,:]-166*y[i+1,:]+170*y[i+2,:]-110*y[i+3,:]+44*y[i+4,:]-10*y[i+5,:]+y[i+6,:])/2;
    else
        delta6=(y[i+3,:]-6*y[i+2,:]+15*y[i+1,:]-20*y[i,:]+15*y[i-1,:]-6*y[i-2,:]+y[i-3,:]);
    end

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

# Finite difference approximation
function Der(y,i,k,X)

    if i==4 # left boundary TEM1
        z = (-27*y[i,k]+58*y[i+1,k]-56*y[i+2,k]+36*y[i+3,k]-13*y[i+4,k]+2*y[i+5,k])/(12*(X[i+1]-X[i]))
    elseif i==5 # left boundary TEM2
        z = (-2*y[i-1,k]-15*y[i,k]+28*y[i+1,k]-16*y[i+2,k]+6*y[i+3,k]-y[i+4,k])/(12*(X[i+1]-X[i]))
    else # central
        z = (-y[i+2,k]+8*y[i+1,k]-8*y[i-1,k]+y[i-2,k])/(12*(X[i+1]-X[i]))
    

    end
    
    return z
end


int(x) = floor(Int, x)


# Test Model RHSs for the bulk equations (3.6.16)


function SFconstraint_beta(beta0,x1,time,funcs)
    
    psi = funcs[1]
    derpsi = funcs[2]
    
    #println("psi(x1)[1]", psi(x1)[1])

    if x1<10^(-15)
        z = 0.0
    elseif abs.(x1 .- 1.0)<10^(-15)
        z = 0.0
    else
        #z = 2.0 .* pi .* (1.0 .- x1) ./ x1 .^3.0 .* (psi(x1) .+ (x1 .- 1) .* x1 .* derpsi(x1)) .^2
        z = 2.0 .* pi .* (1.0 .- x1) ./ x1 .^3.0 .* (psi(x1) .+ (x1 .- 1.0) .* x1 .* derpsi(x1)) .^2.0
        #z = 2.0 .* pi .* (1.0 .- x1) ./ x1 .^3.0 .* ((x1 .- 1.0) .* x1 .* derpsi(x1)[1]) .^2.0
    end

    #z=2.0 .* pi .* (1.0 .- x1) 
    return z
    #return psi(x1)[1]

    #return derpsi(x1)[1]
end

function new_constraint(f0,x1,time,funcs)
    
    z = sin(x1*20+time)

    return z
end


function SFconstraint_m(m0,x1,time,funcs)

    psi = funcs[1]
    derpsi = funcs[2]

    if x1<10^(-15)
        z = 0
    elseif abs.(x1 .- 1.0)<10^(-15)
        z = 2.0 .* pi .* (psi(x1)) .^ 2.0
    else
        #z = 2.0 .* pi .* (x1 .+ 2.0 .* (x1 .- 1.0) .* m0) ./ x1 .^3.0 .* (psiarray[k] .+ (x1 .- 1.0) .* x1 .* derpsi(x1)) .^2.0
        z = 2.0 .* pi .* (x1 .+ 2.0 .* (x1 .- 1.0) .* m0) ./ x1 .^3.0 .* (psi(x1) .+ (x1 .- 1.0) .* x1 .* derpsi(x1)) .^ 2.0
    end
    
    return z
end



function bulkSF(y,i,X)
    
    
    #psi,X

    dy=-1.0/2.0*exp(2.0*y[i,2])*(-(2*(X[i]-1)^3*y[i,3]*(X[i]*((X[i]-1)*Der(y,i,1,X)+X[i]*Der(y,i,2,X))+y[i,1]*(1+2*(X[i]-1)*X[i]*Der(y,i,2,X))))/X[i]^3 - (2*(X[i]-1)^4*(X[i]*((X[i]-1)*Der(y,i,1,X)+X[i]*Der(y,i,2,X))+y[i,1]*(1+2(X[i]-1)*X[i]*Der(y,i,2,X)))*y[i,4])/X[i]^2 - ((X[i]+2*(X[i]-1)*y[i,1])*Der(y,i,4,X))/X[i]) #- dissipation6(y,i,0.01)[4];
    #dy=-1.0/2.0*exp(2.0*y[i,2])*((2.0*exp(2.0*(X[i]-y[i,3]+X[i]*y[i,3])*y[i,2]/X[i])*(X[i]-1.0)^2*(X[i]*((X[i]-1.0)*Der(y,i,1,X)+X[i]*Der(y,i,2,X))+y[i,1]*(1.0+2.0*(X[i]-1.0)*X[i]*Der(y,i,2,X))))/X[i]^2 - ((-1.0+X[i])^3*(X[i]+2.0*(X[i]-1.0)*y[i,1])*y[i,4])/X[i]^2 - ((1.0-X[i])^3*(1.0-2.0*(X[i]-1.0)^2*Der(y,i,1,X))*y[i,4])/X[i] - (2.0*(X[i]-1.0)^4*(X[i]+2.0*(X[i]-1.0)*y[i,1])*Der(y,i,2,X)*y[i,4])/X[i] - (Der(y,i,4,X)) - (2.0*(X[i]-1.0)*y[i,1]*Der(y,i,4,X))/X[i]) #- dissipation6(y,i,0.01)[4];

    return dy
end



function boundarySF(y,X)

    L=length(state_array[:,1])

    dxx=X[L-2]-X[L-3]

    #m even
    y[3,1]=y[5,1];
    y[2,1]=y[6,1];
    y[1,1]=extrapolate_in(y[2,1], y[3,1], y[4,1], y[5,1]);

    y[L-2,1]=y[L-4,1]+dxx/8*pi*(y[L-3,3])^2;
    y[L-1,1]=y[L-5,1]-dxx*pi*(y[L-3,3])^2;
    y[L,1]=extrapolate_out(y[L-4,1], y[L-3,1], y[L-2,1], y[L-1,1]);

    #beta even
    y[3,2]=y[5,2];
    y[2,2]=y[6,2];
    y[1,2]=extrapolate_in(y[2,2], y[3,2], y[4,2], y[5,2]);

    y[L-2,2]=y[L-4,2];
    y[L-1,2]=y[L-5,2];
    y[L,2]=extrapolate_out(y[L-4,2], y[L-3,2], y[L-2,2], y[L-1,2]);
    
    #psi odd
    y[3,3]=-y[5,3];
    y[2,3]=-y[6,3];
    y[1,3]=-y[7,3];

    y[L-2,3]=extrapolate_out(y[L-6,3],y[L-5,3], y[L-4,3], y[L-3,3])
    y[L-1,3]=extrapolate_out(y[L-5,3],y[L-4,3], y[L-3,3], y[L-2,3])
    y[L,3]=extrapolate_out(y[L-4,3], y[L-3,3], y[L-2,3], y[L-1,3]);

    #psi,x even
    y[3,4]=y[5,4];
    y[2,4]=y[6,4];
    y[1,4]=y[7,4];

    y[L-2,4]=extrapolate_out(y[L-6,4],y[L-5,4], y[L-4,4], y[L-3,4])
    y[L-1,4]=extrapolate_out(y[L-5,4],y[L-4,4], y[L-3,4], y[L-2,4])
    y[L,4]=extrapolate_out(y[L-4,4], y[L-3,4], y[L-2,4], y[L-1,4])

    
    return y
end

# Defining the function in the RHS of the evolution equation system

function SF_RHS(data,t,X)
    
    L=length(X)
    dy=zeros((L,length(data[1,:])));

    # update interpolation of psi,x
    spl_derpsi = scipyinterpolate.splrep(X[4:L-3], data[4:L-3,4],k=4)#old
    derpsi_func(x) = scipyinterpolate.splev(x, spl_derpsi)#old
    #derpsi_func = cubic_spline_interpolation(X[4:L-3], data[4:L-3,4],  extrapolation_bc = Line())#new
    
    # rk4wrapper to update psi data
    psi0=0
    SFconstraint_psi(psi0,x) = scipyinterpolate.splev(x, spl_derpsi)#old
    data[4:L-3,3] = rungekutta4(SFconstraint_psi,psi0,X[4:L-3])#old
    #SFconstraint_psi(psi0,x) = derpsi_func(x)#new
    #data[4:L-3,3] = rungekutta4(SFconstraint_psi,psi0,X[4:L-3])#new
    data = ghost(data)
    #global state_array[:,3] = data[:,3]

    # update interpolation of psi
    spl_psi = scipyinterpolate.splrep(X[4:L-3], data[4:L-3,3],k=4)#old
    psi_func(x) = scipyinterpolate.splev(x, spl_psi)#old
    #psi_func = cubic_spline_interpolation(X[4:L-3], data[4:L-3,3],  extrapolation_bc = Line())#new

    funcs = [psi_func derpsi_func]

    # rk4wrapper to update beta data
    beta0=0
    data[4:L-3,2] = rk4wrapper(SFconstraint_beta,beta0,X[4:L-3],t,funcs)
    data = ghost(data)
    #global state_array[:,2] = data[:,2]

    # rk4wrapper to update m data
    m0=0
    data[4:L-3,1]=rk4wrapper(SFconstraint_m,m0,X[4:L-3],t,funcs)
    data = ghost(data)
    #global state_array[:,1] = data[:,1]


 
    #psi another RHS. but do it later




"""
    #identify 1st positive gridpoint X[i]
    origin_i = 4
    for i in 1:L
        if X[i] >= 0
            origin_i = i
            break
        end
    end
    println("\nOrigin of the grid is at t = ", t, " is i = ", origin_i, "\n")
   """

    for i in 4:L-3 #ORI
        if X[i]<10^(-15) #left
            dy[i,4]= +1/2*Der(data,i,4,X)-dissipation6(data,i,0.035)[4];

        elseif X[i] < (1-10^(-15)) #bulk
            dy[i,4]=bulkSF(data,i,X) - dissipation6(data,i,0.035)[4];

        else #right
            dy[i,4]=bulkSF(data,i,X) - dissipation6(data,i,0.035)[4];
        end
    end
    
    
    
    #outer boundary
    #dy[L-3,4]=extrapolate_out(dy[L-7,4], dy[L-6,4], dy[L-5,4], dy[L-4,4])#1.0/2.0*exp(2.0*y[L-3,2])* Der(y,L-3,4)
    #dy[L-2,4]=extrapolate_out(dy[L-6,4], dy[L-5,4], dy[L-4,4], dy[L-3,4])#1.0/2.0*exp(2.0*y[L-3,2])* Der(y,L-2,4)
    #dy[L-1,4]=extrapolate_out(dy[L-5,4], dy[L-4,4], dy[L-3,4], dy[L-2,4])
    #dy[L,4]=extrapolate_out(dy[L-4,4], dy[L-3,4], dy[L-2,4], dy[L-1,4])
  

    return dy

end


function Grid_RHS(y,t,X)
    
    L=length(X)
    dy=zeros((L,length(y[1,:])));

    for i in 4:(L-3)
        if X[i]<10^(-7) #left

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

using ProgressMeter
function timeevolution(state_array,finaltime,dir,dt)

    @showprogress for t in 1:finaltime

        if isnan(state_array[L-3,4])
            print("boom at t=", t*dt)
            explode = true
            timestep = t
            break
        end
        
        #X = update_grid(state_array[:,:],T,t)
        
        X=initX#initX #state_array[:,5]
        X1=X[4:L-3]
    
        #update ghost points
        state_array=boundarySF(state_array,X)
       
        #evolve psi,x
        state_array[:,:] = rungekutta4molstep(SF_RHS,state_array[:,:],T,t,0,X) #evolve psi,x
        #global state_array=ghost(state_array)
    
        #global aux=SF_RHS(state_array[:,:], 0,0,X)
        
        #calculate psi from psi,x
        spl_derpsi = scipyinterpolate.splrep(initX[4:L-3], state_array[4:L-3,4],k=4)
        psi0=0
        SFconstraint_psi(psi0,x) = scipyinterpolate.splev(x, spl_derpsi)
    
        state_array[4:L-3,3] = rungekutta4(SFconstraint_psi,psi0,initX1)
        #global state_array=ghost(state_array);
    
        spl_psi = scipyinterpolate.splrep(initX[4:L-3], state_array[4:L-3,3],k=4)
        psi_func(x) = scipyinterpolate.splev(x, spl_psi)
        derpsi_func(x) = scipyinterpolate.splev(x, spl_derpsi)
    
        funcs = [psi_func derpsi_func]
        
        #evolve beta
        state_array[4:L-3,2]=rk4wrapper(SFconstraint_beta,beta0,X1,T[t+1],funcs)
        #global state_array=ghost(state_array)
        
        #evolve m
        state_array[4:L-3,1]=rk4wrapper(SFconstraint_m,m0,X1,T[t+1],funcs)
        #global state_array=ghost(state_array)
        
        #CSV.write(dir*"/time_step$k.csv", Tables.table(transpose(Matrix(state_array))), writeheader=false)
        
        if t%10==0
            CSV.write(dir*"/time_step$t.csv", Tables.table(state_array), writeheader=false)
        end

        #threshold for apparent black hole formation
        global monitor_ratio = zeros(L)
        
        for i in 1:L
            global monitor_ratio[i] = 2*state_array[i,1]/initX[i]*(1-initX[i])
            if monitor_ratio[i]>1.0
                global criticality = true
                println("Supercritical evolution! At timestep ", t*dt)
                println("i = ", i, " t = ", t, " monitor ratio = ", monitor_ratio[i])
                global timestep = t
            end
        end
        if criticality == true
            break
        end
        

        global timestep = t
    end

    global evol_stats = [criticality A sigma r0 timestep explode]

    CSV.write("/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/DATA/parameters.csv", Tables.table(evol_stats), writeheader=true,header=["criticality", "A", "sigma", "r0", "timestep", "explode"]);

end    