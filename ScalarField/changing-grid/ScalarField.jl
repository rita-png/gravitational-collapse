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

# Calculating dt
function dt_scale(X, m, beta,dx)

    L = length(X)
    g = []
    ori=find_origin(X)
    for i in (ori:(L-4))
        g = vcat(g,exp(2*beta[i])*m[i])
    end
    
    return  maximum(g)
end

"""function update_dt(X, m, beta,dx,ginit)

    monitor_ratio = zeros(L)
    g=ginit

    for i in 1:L
        monitor_ratio[i] = 2*state_array[i,1]/initX[i]*(1-initX[i])
        if monitor_ratio[i]>0.55
            g=dt_scale(X,m,beta,dx)
        end
    end

    return  dx*sqrt(ginit/g)
end
""" #THis needs to be fixed
function find_origin(X)

    #origin_i=Int64
    L=length(X)
    origin_i = 4
    for i in 1:L
        if X[i] >= 0
            origin_i = i
            break
        end
    end
    #println("\nOrigin of the grid is at t = ", t, " is i = ", origin_i, "\n")

    return origin_i::Int64#floor(Int,origin_i)
end

# Updating Grid
function update_grid(data,T,k)
    
    X=data[:,5]
    ori = find_origin(X)
    
    m_func = Spline1D(X[ori:L-3],data[ori:L-3,1],k=4)#new
    beta_func = Spline1D(X[ori:L-3],data[ori:L-3,2],k=4)#new
    psi_func = Spline1D(X[ori:L-3],data[ori:L-3,3],k=4)#new
    derpsi_func = Spline1D(X[ori:L-3],data[ori:L-3,4],k=4)#new

    #evolve grid
    data = rungekutta4molstep(Grid_RHS,data,T,k,data[:,5]) #evolve X here
    data = ghost(data)

    #update X
    X = data[:,5]
    ori = find_origin(X)

    
    #repopulate grid
    if ori>length(X)/2
        println("GRID DUPLICATED!")
        X=doublegrid(X)
        global initX = doublegrid(initX)
        ori = find_origin(X)
    end


    #update m, beta, psi and psi,x data on initial grid
    data[ori:L-3,1]=m_func(initX[ori:L-3])
    data[ori:L-3,2]=beta_func(initX[ori:L-3])
    data[ori:L-3,3]=psi_func(initX[ori:L-3])
    data[ori:L-3,4]=derpsi_func(initX[ori:L-3])

    return X,data[ori:L-3,1:4]

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

function rungekutta4molstep(f,y00,T,w::Int64,X)
    y = y00;
    X=collect(X)
        h = T[w+1]-T[w]
        """#print("\n\nh = ", h, " \n")
        #h = dt; # only for equally spaced grid in time and space, otherwise (T[w+1] - T[w])
        k1 = f(y[:,:], T[w],X)
        #k1=ghost(k1)
        #k1=boundarySF(k1,X)
        k2 = f(y[:,:] .+ k1 .* h/2, T[w] + h/2,X)
        #k2=boundarySF(k2,X)
        k3 = f(y[:,:] .+ k2 .* h/2, T[w] + h/2,X)
        #k3=boundarySF(k3,X)
        k4 = f(y[:,:] .+ k3 .* h, T[w] + h,X)
        #k4=boundarySF(k4,X)
        y[:,:] = y[:,:] .+ (h/6) .* (k1 .+ 2 * k2 .+ 2 * k3 .+ k4)"""
        
        k1 = f(y[:,:], T[w],X)
        k1=ghost(k1)
        k2 = f(y[:,:] .+ k1 .* h, T[w] + h,X)
        k2=ghost(k2)
        y[:,:] = y[:,:] .+ (h/2) .* (k1 .+ k2)
        
    return ghost(y[:,:])
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

function n_rk4wrapper(f,y0,x,u,spl_funcs) # u depicts T array, or M!!
    L = length(x)
    n = length(y0)
    y = zeros(L,n)
    y[1,:] = y0;

    
    """for i in 1:L-1
        h = x[i+1] .- x[i]
        k1 = f(y[i,:], x[i],u,spl_funcs)
        k2 = f(y[i,:] .+ k1 * h/2, x[i] .+ h/2,u,spl_funcs)
        k3 = f(y[i,:] .+ k2 * h/2, x[i] .+ h/2,u,spl_funcs)
        k4 = f(y[i,:] .+ k3 * h, x[i] .+ h,u,spl_funcs)
        y[i+1,:] = y[i,:] .+ (h/6) * (k1 .+ 2*k2 .+ 2*k3 .+ k4)
    end"""
    for i in 1:L-1
        h = x[i+1] .- x[i]
        k1 = f(y[i,:], x[i],u,spl_funcs)
        k2 = f(y[i,:] .+ k1 * h, x[i] .+ h,u,spl_funcs)
        y[i+1,:] = y[i,:] .+ (h/2) * (k1 .+ k2)
    end
    return y[:,:]
end

function print_muninn(io::IO, t, data)
    #@assert length(xs) == length(data[:,1])
    @printf io "\"Time = %.10e\n" t
    for i in 1:length(data[:,1])
        @printf io "% .10e % .10e % .10e % .10e % .10e\n" data[i,5] data[i,1] data[i,2] data[i,3] data[i,4]
    end
    println(io) # insert empty line to indicate end of data set
end


#ghosts

function ghost(y)
    L=length(y[:,1])
    
    #inner boundary extrapolation
    y[3,1:4]=extrapolate_in(y[4,1:4], y[5,1:4], y[6,1:4], y[7,1:4])
    y[2,1:4]=extrapolate_in(y[3,1:4], y[4,1:4], y[5,1:4], y[6,1:4])
    y[1,1:4]=extrapolate_in(y[2,1:4], y[3,1:4], y[4,1:4], y[5,1:4])

    #outer boundary extrapolation
    y[L-2,1:5]=extrapolate_out(y[L-6,1:5], y[L-5,1:5], y[L-4,1:5], y[L-3,1:5])
    y[L-1,1:5]=extrapolate_out(y[L-5,1:5], y[L-4,1:5], y[L-3,1:5], y[L-2,1:5])
    y[L,1:5]=extrapolate_out(y[L-4,1:5], y[L-3,1:5], y[L-2,1:5], y[L-1,1:5])
   

    return y
end


#6th order dissipation, added to 4th order original scheme
function dissipation6(y,i,eps)

    ori = find_origin(y[:,5])
    if i==ori
        delta6= (19*y[i,:]-142*y[i+1,:]+464*y[i+2,:]-866*y[i+3,:]+1010*y[i+4,:]-754*y[i+5,:]+352*y[i+6,:]-94*y[i+7,:]+11*y[i+8,:])/2;
    elseif i==ori+1
        delta6= (11*y[i-1,:]-80*y[i,:]+254*y[i+1,:]-460*y[i+2,:]+520*y[i+3,:]-376*y[i+4,:]+170*y[i+5,:]-44*y[i+6,:]+5*y[i+7,:])/2;
    elseif i==ori+2
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


#2nd order  dissipation, added to 1st order scheme
function dissipation2(y,i)
        delta2=(y[i+1,:]-2*y[i,:]+y[i-1,:]);
    return (-1)^1*epsilon*1/(dx)*delta2
end


# Finite difference approximation
function Der(y,i,k,X)
    ori = find_origin(X)
    if i==ori # left boundary TEM1
        z = (-27*y[i,k]+58*y[i+1,k]-56*y[i+2,k]+36*y[i+3,k]-13*y[i+4,k]+2*y[i+5,k])/(12*(X[i+1]-X[i]))
    elseif i==ori+1 # left boundary TEM2
        z = (-2*y[i-1,k]-15*y[i,k]+28*y[i+1,k]-16*y[i+2,k]+6*y[i+3,k]-y[i+4,k])/(12*(X[i+1]-X[i]))
    else # central
        z = (-y[i+2,k]+8*y[i+1,k]-8*y[i-1,k]+y[i-2,k])/(12*(X[i+1]-X[i]))
    end
    
    return z
end

function Der(y,i,X) #dx tilde/dx
    ori = find_origin(y)
    if i==ori # left boundary TEM1
        z = (-27*y[i]+58*y[i+1]-56*y[i+2]+36*y[i+3]-13*y[i+4]+2*y[i+5])/(12*(X[i+1]-X[i]))
    elseif i==ori+1 # left boundary TEM2
        z = (-2*y[i-1]-15*y[i]+28*y[i+1]-16*y[i+2]+6*y[i+3]-y[i+4])/(12*(X[i+1]-X[i]))
    else # central
        z = (-y[i+2]+8*y[i+1]-8*y[i-1]+y[i-2])/(12*(X[i+1]-X[i]))
    

    end
    
    return z
end


int(x) = floor(Int, x)


function bulkSF(y,i,X)
    
    
    #psi,X

    dy=-1.0/2.0*exp(2.0*y[i,2])*(-(2*(X[i]-1)^3*y[i,3]*(X[i]*((X[i]-1)*Der(y,i,1,X)+X[i]*Der(y,i,2,X))+y[i,1]*(1+2*(X[i]-1)*X[i]*Der(y,i,2,X))))/X[i]^3 - (2*(X[i]-1)^4*(X[i]*((X[i]-1)*Der(y,i,1,X)+X[i]*Der(y,i,2,X))+y[i,1]*(1+2(X[i]-1)*X[i]*Der(y,i,2,X)))*y[i,4])/X[i]^2 - ((X[i]+2*(X[i]-1)*y[i,1])*Der(y,i,4,X))/X[i]) #- dissipation6(y,i,0.01)[4];
    
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


function RHS(y0,x1,time,func)
    
    z=zeros(length(y0))
    derpsi = func[1]
    der_grid = func[2]

    z[3] = derpsi(x1)

    
    if derpsi(x1)>10
        println("  ")
        println(derpsi(x1))
        println("x1 isssssssss ",x1)
    end#IF X1<0, extrapolate_in

    if abs.(x1)<10^(-15) #origin
        z[1] = 0.0;
        z[2] = 0.0;
    elseif abs.(x1 .- 1.0)<10^(-15) #right
        z[1] = 2.0 .* pi .* (y0[3]) .^ 2.0 ./ der_grid(x1)
        z[2] = 0.0
    else #right
        z[1] = 2.0 .* pi .* (x1 .+ 2.0 .* (x1 .- 1.0) .* y0[1]) ./ x1 .^3.0 .* (y0[3] .+ (x1 .- 1.0) .* x1 .* derpsi(x1)) .^ 2.0 ./ der_grid(x1);
        z[2] = 2.0 .* pi .* (1.0 .- x1) ./ x1 .^3.0 .* (y0[3] .+ (x1 .- 1.0) .* x1 .* derpsi(x1)) .^2.0 ./ der_grid(x1);
    end


    return z[:]
end

# freely falling points 
function null_ingoing(y,i,X)

    z = -1.0/2.0*(1-X[i])^2.0*exp(2*y[i,2])*(1-2*y[i,1]*(1-X[i])/X[i])

    return z
end


function leftboundary(data,funcs)
    xtilde_func=Spline1D(initX1,data[4:L-3,5]);
    println("xtilde_func(0) ",xtilde_func(0))
    
    ori=find_origin(data[:,5])
    auxX=vcat(xtilde_func(0),data[ori,5]);

    println("auxX",auxX)
    
    y0=[0 0 0]
    auxdata=zeros(2,3)
    auxdata[:,1:3] = n_rk4wrapper(RHS,y0,auxX,0,funcs) #problem, these funcs are not defined for the negative x i need here
    println("auxdata[:,1:3]", auxdata[:,1:3])
    return auxdata[2,1:3]
end


# Defining the function in the RHS of the evolution equation system

function SF_RHS(data,t,X)
    
    L=length(X)
    dy=zeros((L,length(data[1,:])));
    ori = find_origin(X)

    # update interpolation of psi,x
    derpsi_func = Spline1D(X[ori:L-3],data[ori:L-3,4],k=4, bc="extrapolate")#new
    dergrid_func=der_grid(X)
    funcs=[derpsi_func dergrid_func]

    ##
    """xtilde_func=Spline1D(initX1,data[4:L-3,5]);
    xtilde_func(0)
    
    auxX=vcat(xtilde_func(0),X[ori:L-3]);"""
    """println("  ")
    println(auxX)"""
    ##

    #println(leftboundary(data,funcs))
    # update m, beta and psi data
    y0=[0 0 0]
    #y0=[0.0000001 0.0000001 0.0000001]#leftboundary(data,funcs)
    data[ori:L-3,1:3] = n_rk4wrapper(RHS,y0,X[ori:L-3],t,funcs)


    for i in ori:L-3
        if X[i]<10^(-15) #left
            dy[i,4]= +1/2*Der(data,i,4,X) - null_ingoing(data,i,X)*Der(data,i,4,initX)/(Der(X,i,initX)) - dissipation6(data,i,0.03)[4]; #This RHS is now evaluated at x tilde, not x

        elseif X[i] < (1-10^(-15)) #bulk
            dy[i,4]=bulkSF(data,i,X) - null_ingoing(data,i,X)*Der(data,i,4,initX)/(Der(X,i,initX)) - dissipation6(data,i,0.03)[4];

        else #right
            dy[i,4]=bulkSF(data,i,X) - null_ingoing(data,i,X)*Der(data,i,4,initX)/(Der(X,i,initX)) - dissipation6(data,i,0.03)[4];
        end
    end

    return dy

end


function Grid_RHS(y,t,X)
    
    L=length(X)
    dy=zeros((L,length(y[1,:])));

    for i in 4:(L-3)
        if X[i]<10^(-15) && X[i]>-10^(-15) #around origin

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

#calculate dx~/dx
function der_grid(X)

    aux=[]
    ori=find_origin(X)

    for i in ori:L-3
        aux = vcat(aux,Der(X,i,initX))
    end
    
    dergrid_func = Spline1D(X[ori:L-3],aux,k=4, bc="extrapolate")

    return dergrid_func
end
    


using Term.Progress
function timeevolution(state_array,finaltime,dir,dt,run)

    t=0.0
    T_array = [0.0]
    T_interp = [0.0]
    iter = 0

    while t<=finaltime#@TRACK

        iter = iter + 1

        #update time increment
        #global dt = update_dt(state_array[:,5],state_array[:,1],state_array[:,2],dx,ginit)
        t = t + round(dt,digits=5)
        """if iter%10==0
            println("iteration ", iter, " dt is ", dt, ", time of iteration is ", t)
        end"""
        println("iteration ", iter, " dt is ", dt, ", time of iteration is ", t)

        T_array = vcat(T_array,t)

        #update grid
        X,aux = update_grid(state_array[:,:],T_array,iter)
        ori=find_origin(X)
        state_array[ori:L-3,1:4]=aux[:,1:4]
        state_array[:,5]=X

        
        X1=X[ori:L-3]
        
        #evolve psi,x
        state_array[:,:] = rungekutta4molstep(SF_RHS,state_array[:,:],T_array,iter,X) #evolve psi,x using data on initX grid
        #global state_array=ghost(state_array)
    
        # update interpolation of psi,x
        derpsi_func = Spline1D(X[ori:L-3],state_array[ori:L-3,4],k=4,bc="extrapolate")#new #OLAOLA
        dergrid_func = der_grid(X)
        funcs=[derpsi_func dergrid_func]
        
        ##
        #leftboundary(state_array,funcs)
        """xtilde_func=Spline1D(initX1,state_array[4:L-3,5]);
        xtilde_func(0)

        auxX=vcat(xtilde_func(0),X[ori:L-3]);"""
        
        ##

        # update m, beta and psi data
        y0=[0 0 0]
        #y0=leftboundary(state_array,funcs) THIS SHOULD BE UNCOMMENTED BUT IT CRASHES WHEN I UNCOMMENT IT
        println("timestep ", leftboundary(state_array,funcs))
        state_array[ori:L-3,1:3] = n_rk4wrapper(RHS,y0,X[ori:L-3],t,funcs) #OLAOLA
        

        CSV.write(dir*"/res$res/time_step$iter.csv", Tables.table(state_array), writeheader=false)
        run=int(run)
        """if iter%10==0
            #CSV.write(dir*"/run$run/time_step$iter.csv", Tables.table(state_array), writeheader=false)
            CSV.write(dir*"/res$res/time_step$iter.csv", Tables.table(state_array), writeheader=false)
            T_interp = vcat(T_interp,t)

            #write muninn
            open(dir*"/res$res/data.txt", "a") do file
                print_muninn(file, t, state_array[:,1:5])
            end
        end"""
        
        

        #threshold for apparent black hole formation
        global monitor_ratio = zeros(L)
        
        for i in 4:L-3
            global monitor_ratio[i] = 2*state_array[i,1]/X[i]*(1-X[i])
            if monitor_ratio[i]>1.0
                global criticality = true
                println("Supercritical evolution! At time ", t)
                println("Gridpoint = ", i, " t = ", t, " monitor ratio = ", monitor_ratio[i])
                global time = t
            end
        end
        
        if criticality == true
            break
        end
        
        if isnan(state_array[L-3,4])
            global explode = true
            println("boom at time=", t)
            global time = t
            break
        end

        global time = t
        
    end
    
    global evol_stats = [criticality A sigma r0 time explode run]

    return evol_stats,T_interp
   
end    