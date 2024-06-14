# Definition of gaussian initial data functions

function init_derpsi(r,r0,sigma,A) #psi,r
    n=length(r);
    if compactified==false
        if n==1
            z= A*exp(-((r - r0)^2/sigma^2)) - (2*A*exp(-((r - r0)^2/sigma^2))*r*(r - r0))/sigma^2#A * (3 * exp(-(r-r0)^2/sigma^2) * r^2 - 2 * exp(-(r-r0)^2/sigma^2) * (r-r0)*r^3/sigma^2)#exp(-((x/(1-x)-r0)/sigma)^2) * (3 * x^2 / (1-x) ^4 - (x/(1-x))^3 * (2*(x-r0*(-x+1)))/(sigma^2*(1-x)^3))
        else
            z=zeros(n);
            for i in 1:n
                rr = r[i]
                z[i] = A*exp(-((rr - r0)^2/sigma^2)) - (2*A*exp(-((rr - r0)^2/sigma^2))*rr*(rr - r0))/sigma^2#A * (3 * exp(-(rr-r0)^2/sigma^2) * rr^2 - 2 * exp(-(rr-r0)^2/sigma^2) * (rr-r0)*rr^3/sigma^2)
            end
        end
    else # inputted argument r is actually an x
        if loggrid==false
            if n==1
                x=r
                rr=x/(1-x)
                z= A*exp(-((rr - r0)^2/sigma^2)) - (2*A*exp(-((rr - r0)^2/sigma^2))*rr*(rr - r0))/sigma^2#A * (3 * exp(-(rr-r0)^2/sigma^2) * rr^2 - 2 * exp(-(rr-r0)^2/sigma^2) * (rr-r0)*rr^3/sigma^2)
            else
                z=zeros(n);
                for i in 1:n

                    ## psi,r
                    x=r[i]
                    rr = x/(1-x)
                    z[i] = A*exp(-((rr - r0)^2/sigma^2)) - (2*A*exp(-((rr - r0)^2/sigma^2))*rr*(rr - r0))/sigma^2#A * (3 * exp(-(rr-r0)^2/sigma^2) * rr^2 - 2 * exp(-(rr-r0)^2/sigma^2) * (rr-r0)*rr^3/sigma^2)

                    
                end
                z[n] = 0
            end
        else
            z=0
        end
    end

    return z
end

function init_derxi(r,r0,sigma,A)
    n=length(r);
    if compactified==false
        if n==1
            z= 0#?
        else
            z=zeros(n);
            for i in 1:n
                rr = r[i]
                z[i] = 0#?
            end
        end
    else # inputted argument r is actually an x
        
        if n==1
            x=r
            rr=x/(1-x)
            z= -((2*A*exp(-((rr+r0)^2/sigma^2))*rr*((1+exp((4*rr*r0)/sigma^2))*rr^2-(-1+exp((4*rr*r0)/sigma^2))*rr*r0-(1+exp((4*rr*r0)/sigma^2))*sigma^2))/sigma^2)
        else
            z=zeros(n);
            for i in 1:n

                ## psi,r
                x=r[i]
                rr = x/(1-x)
                z[i] = -((2*A*exp(-((rr+r0)^2/sigma^2))*rr*((1+exp((4*rr*r0)/sigma^2))*rr^2-(-1+exp((4*rr*r0)/sigma^2))*rr*r0-(1+exp((4*rr*r0)/sigma^2))*sigma^2))/sigma^2)

                    
            end
       
        end
    end

    return z
end
    

# Extrapolation

function extrapolate_out(y0,y1,y2,y3)
    """if length(y0)==1
        println("len is 1")
        return -y0 + 4*y1 - 6*y2 + 4*y3
    else
        println("len is not 1")
        z=zeros(length(y0))
        for i in 1:length(y0)
            z[i]=y0[i] + 4*y1[i] - 6*y2[i] + 4*y3[i]
        end
        ?
        return z
    end"""
    return -y0 + 4*y1 - 6*y2 + 4*y3
end

function extrapolate_in(y0,y1,y2,y3)
    return -y3 + 4*y2 - 6*y1 + 4*y0
end



# Calculating dt
function speed(X, m, beta)

    L = length(X)
    
    ori=find_origin(X)

    g = zeros(int(L-5-ori))
    #g=abs.((1.0 .- initX[ori+1:L-4]) .^ 3.0 .* exp.(2 .* state_array[ori+1:L-4,2]) .* (2 .* state_array[ori+1:L-4,1] .- initX[ori+1:L-4] ./ (1 .- initX[ori+1:L-4])) ./ (2 .* initX[ori+1:L-4]))
    
    if compactified==true
        g=abs.((1.0 .- initX[ori+1:L-4]) .^ 3.0 .* exp.(2 .* state_array[ori+1:L-4,2]) .* (2 .* state_array[ori+1:L-4,1] .- initX[ori+1:L-4] ./ (1 .- initX[ori+1:L-4])) ./ (2 .* initX[ori+1:L-4]))
    else
        g=abs.(exp.(2 .* state_array[ori+1:L-4,2]) .* (2 .* state_array[ori+1:L-4,1] .- initX[ori+1:L-4]) ./ (2 .* initX[ori+1:L-4]))
    end
    z=maximum(g)
    if isnan(z)
        println("Error: Speed is NaN!")
    end
    return z
end


function update_dt(X, m, beta,dt,ginit)

    monitor_ratio = zeros(L)
    
    g=speed(X,m,beta)

    if loggrid==false
        dx=X[5]-X[4]
    else
        aux=zeros(L-7)
        for i in 1:L-7
            aux[i]=initX1[i+1]-initX1[i]
        end
        
        dx=minimum(aux)
    end
    return dx/g*0.1
   

end

function find_origin(X)

    origin_i=Int64
    L=length(X)
    
    for i in 1:L
        if X[i] >= 0
            origin_i = i
            break
        end
    end
   
    return origin_i
end

# Updating Grid
function update_grid(data,T,k)
    
    dx = data[2,5]-data[1,5]

    #evolve grid
    data = twod_rungekutta4molstep(Grid_RHS,data,T,k,data[:,5]) #evolve X
    #data=ghost(data)
    
    X = data[:,5]

    

    return X

end
    
#Building initial data with a Runge-Kutta integrator for the constraint

function rungekutta4(f,y0,T)
    n = length(T)
    y = zeros(n)
    y[1] = y0;
    for i in 1:n-1
        h = T[i+1] .- T[i]
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
        #print("\n\nh = ", h, " \n")
        #h = dt; # only for equally spaced grid in time and space, otherwise (T[w+1] - T[w])
        k1 = f(y[:,:], T[w],X)
        k1=ghost(k1)
        k2 = f(y[:,:] .+ k1 .* h/2, T[w] + h/2,X)
        k2=ghost(k2)
        k3 = f(y[:,:] .+ k2 .* h/2, T[w] + h/2,X)
        k3=ghost(k3)
        k4 = f(y[:,:] .+ k3 .* h, T[w] + h,X)
        k4=ghost(k4)
        y[:,:] = y[:,:] .+ (h/6) .* (k1 .+ 2 * k2 .+ 2 * k3 .+ k4)
        
    return ghost(y[:,:])
end

function twod_rungekutta4molstep(f,y00,T,w::Int64,X)
    y = y00;
    X=collect(X)
        h = T[w+1]-T[w]
        #print("\n\nh = ", h, " \n")
        #h = dt; # only for equally spaced grid in time and space, otherwise (T[w+1] - T[w])
        k1 = f(y[:,:], T[w],X)
        k1=ghost(k1)
        k2 = f(y[:,:] .+ k1 .* h, T[w] + h,X)
        k2=ghost(k2)
        y[:,:] = y[:,:] .+ (h/2) .* (k1 .+ k2)
        
    return ghost(y[:,:])
end


function rk4wrapper(f,y0,x,u,spl_funcs,data) # u depicts T array, or M!!
    n = length(x)
    y = zeros(n)
    y[1] = y0;
    for i in 1:n-1
        h = x[i+1] .- x[i]
        k1 = f(y[i], x[i],u,spl_funcs,i,data)
        k2 = f(y[i] .+ k1 * h/2, x[i] .+ h/2,u,spl_funcs,i,data)
        k3 = f(y[i] .+ k2 * h/2, x[i] .+ h/2,u,spl_funcs,i,data)
        k4 = f(y[i] .+ k3 * h, x[i] .+ h,u,spl_funcs,i,data)
        y[i+1] = y[i] .+ (h/6) * (k1 .+ 2*k2 .+ 2*k3 .+ k4)
    end
    return y
end

function n_rk4wrapper(f,y0,x,u,spl_funcs,data) # u depicts T array, or M!!
    L = length(x)
    n = length(y0)
    y = zeros(L,n)
    y[1,:] = y0;
    for i in 1:L-1
        h = x[i+1] .- x[i]
        k1 = f(y[i,:], x[i],u,spl_funcs,i,data)
        k2 = f(y[i,:] .+ k1 * h/2, x[i] .+ h/2,u,spl_funcs,i,data)
        k3 = f(y[i,:] .+ k2 * h/2, x[i] .+ h/2,u,spl_funcs,i,data)
        k4 = f(y[i,:] .+ k3 * h, x[i] .+ h,u,spl_funcs,i,data)
        y[i+1,:] = y[i,:] .+ (h/6) * (k1 .+ 2*k2 .+ 2*k3 .+ k4)
    end
    return y[:,:]
end

function twod_n_rk4wrapper(f,y0,x,u,spl_funcs,data) # u depicts T array, or M!!
    L = length(x)
    n = length(y0)
    y = zeros(L,n)
    y[1,:] = y0;
    for i in 1:L-1
        h = x[i+1] .- x[i]
        k1 = f(y[i,:], x[i],u,spl_funcs,i,data)
        k2 = f(y[i,:] .+ k1 * h, x[i] .+ h,u,spl_funcs,i,data)
        y[i+1,:] = y[i,:] .+ (h/2) * (k1 .+ k2)
    end
    return y[:,:]
end


using Printf
function print_muninn(files, t, data, res, mode, grid)
    #mode is "a" for append or "w" for write
    j=1
    #LL=length(state_array[1,:])
    if bisection==false
        for fl in files #normal run
            
            open(dir*"/muninnDATA/res$res/$fl.txt", mode) do file
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


function print_monitorratio(filename, t, data, mode, grid)
    #mode is "a" for append or "w" for write
    j=1
    #LL=length(state_array[1,:])
    
    if loggrid==true
        auxdir= dir*"/bisectionsearch/muninnDATA/uneven"
    else
        auxdir= dir*"/bisectionsearch/muninnDATA/even"
    end
    
    
        
    open(auxdir*"/run$run/$filename.txt", mode) do file
        @printf file "\"Time = %.10e\n" t
        for i in 1:length(data)
            @printf file "% .10e % .10e\n" grid[i] data[i]
        end
        println(file) # insert empty line to indicate end of data set
        end
        
end



# 0 dimension output, save every variable at the ori and scri+
function zero_print_muninn(files, t, data, res, mode)
    #mode is "a" for append or "w" for write
    j=1
    
    if bisection==false
        for fl in files #normal run
            
            open(dir*"/muninnDATA/res$res/$fl.txt", mode) do file                
                @printf file "% .10e % .10e % .10e\n" t data[4,j] data[L-3,j]
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
                @printf file "% .10e % .10e % .10e\n" t data[4,j] data[L-3,j]
            end
            j=j+1
        end
    end
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

function ghostscri(y)
    L=length(y[:,1])
    nvars=length(y[1,:])
    
    #outer boundary extrapolation
    y[L-3,1:nvars]=extrapolate_out(y[L-7,1:nvars], y[L-6,1:nvars], y[L-5,1:nvars], y[L-4,1:nvars])
   
    return y
end

function parity(xchi)

    L=length(xchi)
    xchi[1]=-xchi[7]
    xchi[2]=-xchi[6]
    xchi[3]=-xchi[5]
    xchi[4]=0

    xchi[L-2]=extrapolate_out(xchi[L-6], xchi[L-5], xchi[L-4], xchi[L-3])
    xchi[L-1]=extrapolate_out(xchi[L-5], xchi[L-4], xchi[L-3], xchi[L-2])
    xchi[L]=extrapolate_out(xchi[L-4], xchi[L-3], xchi[L-2], xchi[L-1])

    return xchi
end

function secondparity(derxchi)

    L=length(derxchi)
    derxchi[1]=derxchi[7]
    derxchi[2]=derxchi[6]
    derxchi[3]=derxchi[5]

    derxchi[L-2]=extrapolate_out(derxchi[L-6], derxchi[L-5], derxchi[L-4], derxchi[L-3])
    derxchi[L-1]=extrapolate_out(derxchi[L-5], derxchi[L-4], derxchi[L-3], derxchi[L-2])
    derxchi[L]=extrapolate_out(derxchi[L-4], derxchi[L-3], derxchi[L-2], derxchi[L-1])

    return derxchi
end

function dissipation(y,i,eps,var)

    if twod==true
        return dissipation4(y,i,eps,var)
    else
        return dissipation6(y,i,eps,var) 
    end
end

#6th order dissipation, added to 4th order original scheme
function dissipation6(y,i,eps,var)
    if i==4
        delta6= (19*y[i,:]-142*y[i+1,:]+464*y[i+2,:]-866*y[i+3,:]+1010*y[i+4,:]-754*y[i+5,:]+352*y[i+6,:]-94*y[i+7,:]+11*y[i+8,:])/2;
    elseif i==5
        delta6= (11*y[i-1,:]-80*y[i,:]+254*y[i+1,:]-460*y[i+2,:]+520*y[i+3,:]-376*y[i+4,:]+170*y[i+5,:]-44*y[i+6,:]+5*y[i+7,:])/2;
    elseif i==6
        delta6= (5*y[i-2,:]-34*y[i-1,:]+100*y[i,:]-166*y[i+1,:]+170*y[i+2,:]-110*y[i+3,:]+44*y[i+4,:]-10*y[i+5,:]+y[i+6,:])/2;
    elseif i==L-3
        delta6= (19*y[i,:]-142*y[i-1,:]+464*y[i-2,:]-866*y[i-3,:]+1010*y[i-4,:]-754*y[i-5,:]+352*y[i-6,:]-94*y[i-7,:]+11*y[i-8,:])/2;
    elseif i==L-4
        delta6= (11*y[i+1,:]-80*y[i,:]+254*y[i-1,:]-460*y[i-2,:]+520*y[i-3,:]-376*y[i-4,:]+170*y[i-5,:]-44*y[i-6,:]+5*y[i-7,:])/2;
    elseif i==L-5
        delta6= (5*y[i+2,:]-34*y[i+1,:]+100*y[i,:]-166*y[i-1,:]+170*y[i-2,:]-110*y[i-3,:]+44*y[i-4,:]-10*y[i-5,:]+y[i-6,:])/2;
    else
        delta6=(y[i+3,:]-6*y[i+2,:]+15*y[i+1,:]-20*y[i,:]+15*y[i-1,:]-6*y[i-2,:]+y[i-3,:]);
    end

return (-1)^3*eps*1/(dx)*delta6
end


#4th order  dissipation, added to 2nd order original scheme
function dissipation4(y,i,eps,var)#0.02

    if i==4
        delta4=(-13/6*y[i+5,:]+71/6*y[i+4,:]-77/3*y[i+3,:]+83/3*y[i+2,:]-89/6*y[i+1,:]+19/6*y[i,:])
    elseif i==5
        delta4=(-7/6*y[i+4,:]+41/6*y[i+3,:]-47/3*y[i+2,:]+53/3*y[i+1,:]-59/6*y[i,:]+13/6*y[i-1,:])
    elseif i==L-3
        delta4=-(-13/6*y[i-5,:]+71/6*y[i-4,:]-77/3*y[i-3,:]+83/3*y[i-2,:]-89/6*y[i-1,:]+19/6*y[i,:])
    elseif i==L-4
        delta4=-(-7/6*y[i-4,:]+41/6*y[i-3,:]-47/3*y[i-2,:]+53/3*y[i-1,:]-59/6*y[i,:]+13/6*y[i+1,:])
    else
        delta4=(y[i+2,:]-4*y[i+1,:]+6*y[i,:]-4*y[i-1,:]+y[i-2,:]);
    end
    
    return (-1)^2*eps*1/(dx)*delta4
end

# Finite difference approximation
function Der(y,i,k,X)

    if twod == false
        if i==4 # left boundary LOP2, TEM
            z = (-27*y[i,k]+58*y[i+1,k]-56*y[i+2,k]+36*y[i+3,k]-13*y[i+4,k]+2*y[i+5,k])/(12*(X[i+1]-X[i]))
        elseif i==5 # left boundary LOP1, TEM
            z = (-2*y[i-1,k]-15*y[i,k]+28*y[i+1,k]-16*y[i+2,k]+6*y[i+3,k]-y[i+4,k])/(12*(X[i+1]-X[i]))
        elseif i==L-3 # right boundary TEM
            z = (-27*y[i,k]+58*y[i-1,k]-56*y[i-2,k]+36*y[i-3,k]-13*y[i-4,k]+2*y[i-5,k])/(12*(X[i]-X[i-1])) #fixed this *-1
        elseif i==L-4 # right boundary TEM
            z = (-2*y[i+1,k]-15*y[i,k]+28*y[i-1,k]-16*y[i-2,k]+6*y[i-3,k]-y[i-4,k])/(12*(X[i+1]-X[i]))
        else # central
            z = (-y[i+2,k]+8*y[i+1,k]-8*y[i-1,k]+y[i-2,k])/(12*(X[i+1]-X[i]))
        end
    else
        if i==4 # left boundary LOP1, TEM
            z = (y[i+3,k]-4*y[i+2,k]+7*y[i+1,k]-4*y[i,k])/(2*(X[i+1]-X[i]))
        elseif i==L-3
            z = (-y[i-3,k]+4*y[i-2,k]-7*y[i-1,k]+4*y[i,k])/(2*(X[i]-X[i-1]))
        else
            z = (y[i+1,k]-y[i-1,k])/(2*(X[i+1]-X[i]))
        end
    end

    return z
end


# Finite difference approximation
function DDer(y,i,k,X)

    if i==4 # left boundary LOP2, TEM
        z = (15/4*y[i,k]-77/6*y[i+1,k]+107/6*y[i+2,k]-13*y[i+3,k]+61/12*y[i+4,k]-5/6*y[i+5,k])/((X[i+1]-X[i]))
    elseif i==5 # left boundary LOP1, TEM
        z = (-y[i+5,k]+7*y[i+4,k]-21*y[i+3,k]+34*y[i+2,k]-19*y[i+1,k]-9*y[i,k]+9*y[i-1,k])/(12*(X[i+1]-X[i]))
    else # central
        z = (-y[i+2,k]+16*y[i+1,k]-30*y[i,k]+16*y[i-1,k]-y[i-2,k])/(12*(X[i+1]-X[i]))
    

    end
    
    return z
end
function DDer_array(y,k,X)#4th order accurate 2nd der

    z=zeros(L)
    for i in 4:L-3
        if i==4 # left boundary LOP2, TEM
            z[i] = (15/4*y[i,k]-77/6*y[i+1,k]+107/6*y[i+2,k]-13*y[i+3,k]+61/12*y[i+4,k]-5/6*y[i+5,k])/((X[i+1]-X[i]))
        elseif i==5 # left boundary LOP1, TEM
            z[i] = (-y[i+5,k]+7*y[i+4,k]-21*y[i+3,k]+34*y[i+2,k]-19*y[i+1,k]-9*y[i,k]+9*y[i-1,k])/(12*(X[i+1]-X[i]))
        elseif i==L-3
            z[i] = (15/4*y[i,k]-77/6*y[i-1,k]+107/6*y[i-2,k]-13*y[i-3,k]+61/12*y[i-4,k]-5/6*y[i-5,k])/((X[i+1]-X[i]))
        elseif i==L-4
            z[i] = (-y[i-5,k]+7*y[i-4,k]-21*y[i-3,k]+34*y[i-2,k]-19*y[i-1,k]-9*y[i,k]+9*y[i+1,k])/(12*(X[i+1]-X[i]))
        else # central
            z[i] = (-y[i+2,k]+16*y[i+1,k]-30*y[i,k]+16*y[i-1,k]-y[i-2,k])/(12*(X[i+1]-X[i]))
        end

    end
    
    return z
end

function Der_arrayLOP(y,k,X)#returns darray/dr or darray/dx if non compactified

    z=zeros(length(y[:,k]))
 

    for i in 4:L-3
        if i==4 # left boundary LOP1, TEM
            z[i] = (y[i+3,k]-4*y[i+2,k]+7*y[i+1,k]-4*y[i,k])/(2*(X[i+1]-X[i]))
        elseif i==L-3
            z[i] = (-y[i-3,k]+4*y[i-2,k]-7*y[i-1,k]+4*y[i,k])/(2*(X[i]-X[i-1]))
        else
            z[i] = (y[i+1,k]-y[i-1,k])/(2*(X[i+1]-X[i]))
        end
        
    end
    
    return z
end



int(x) = floor(Int, x)

function chebyshev(N)

    X=zeros(N)
    
    for i in 1:N
        #X[i]=1/2+1/2*cos((2*i-1)*pi/(2*N))
        if i==1
            X[i]=0.0
        else
            X[i]=1/2+1/2*cos((2*i-1)*pi/(2*N))
        end
    end

    return sort(X)
end



#solving for xi,x
function derxi_evol(y,i,X)
    
    #psi,r
    if compactified == false
        r=X[i]
        dy=0.0#?

    else
        x=X[i]
        #this is gooddy=((1/(2*x^2))*exp(2*y[i,2])*(-1+x)^2*(1+y[i,3]-(1+y[i,3])^3+2*(-1+x)^2*y[i,4]*(x*((-1+x)*Der(y,i,1,X)+x*Der(y,i,2,initX))+y[i,1]*(1+2*(-1+x)*x*Der(y,i,2,initX)))+(x*(x+2*(-1+x)*y[i,1])*Der(y,i,4,initX))/(-1+x)^2))

        #same as good but y[i,4]->Der(y,i,3,X)
        dy=((1/(2*x^2))*exp(2*y[i,2])*(-1+x)^2*(1+y[i,3]-(1+y[i,3])^3+2*(-1+x)^2*Der(y,i,3,X)*(x*((-1+x)*Der(y,i,1,X)+x*Der(y,i,2,initX))+y[i,1]*(1+2*(-1+x)*x*Der(y,i,2,initX)))+(x*(x+2*(-1+x)*y[i,1])*((x-1)^2*Der(y,i,4,X)))/(-1+x)^2))

    end
    return dy
end

function hessian_control(data, t)

    y=data[:,4]
    x=initX

    dx=x[6]-x[5]

    control=100000
    result=false

    hess=zeros(length(data[:,4]))
    for i in 5:L-4 #excluding x=0 and x=1
        hess[i]=abs((y[i+1]-2*y[i]+y[i-1])/dx^2)
    end
    if maximum(hess)>control
        result=true
        println("\n\nHessian is really big!\n\n", maximum(hess), " time is ", t)
    end   
    


    return result

end


function boundarySF(y,X)

    L=length(state_array[:,1])
    
    #m even
    y[L-2,1]=y[L-4,1]+dxx/8*pi*(y[L-3,3])^2;
    y[L-1,1]=y[L-5,1]-dxx*pi*(y[L-3,3])^2;
    y[L,1]=extrapolate_out(y[L-4,1], y[L-3,1], y[L-2,1], y[L-1,1]);

    #beta even
    y[L-2,2]=y[L-4,2];
    y[L-1,2]=y[L-5,2];
    y[L,2]=extrapolate_out(y[L-4,2], y[L-3,2], y[L-2,2], y[L-1,2]);
    
    #psi odd
    y[L-2,3]=extrapolate_out(y[L-6,3],y[L-5,3], y[L-4,3], y[L-3,3])
    y[L-1,3]=extrapolate_out(y[L-5,3],y[L-4,3], y[L-3,3], y[L-2,3])
    y[L,3]=extrapolate_out(y[L-4,3], y[L-3,3], y[L-2,3], y[L-1,3]);

    #psi,x even
    y[L-2,4]=extrapolate_out(y[L-6,4],y[L-5,4], y[L-4,4], y[L-3,4])
    y[L-1,4]=extrapolate_out(y[L-5,4],y[L-4,4], y[L-3,4], y[L-2,4])
    y[L,4]=extrapolate_out(y[L-4,4], y[L-3,4], y[L-2,4], y[L-1,4])

    
    return y
end



function RHS(y0,x1,time,func,i,data)
    
    #state array is [m beta xi derxi]
    z=zeros(length(y0))
    
    derxi = func[1]

    #coords
    r=0
    if compactified==false
        r=x1
    else
        x=x1
        r=x/(1-x)
    end

    #xi
    if compactified==true
        
        z[3]=derxi(x1)/(1-x1)^2 #xi,x=xi,r*dr/dx
        if abs.(x1-1)<10^(-15)
            z[3]=0
        end
    else
        z[3]=0#?
    end
    
    
    #beta
    if compactified==true
        if abs.(x1)<10^(-15)
            z[2]=0
        else
            z[2]=-((4*pi*(-1+x1)^3*(z[3])^2)/x1)
        end
    else
        z[2]=0#?
    end

    #m
    if compactified==true
        if abs.(x1)<10^(-15)
            z[2]=0
        else
            x=x1
            z[1]=(2*pi*(y0[3]^2*(2+y0[3])^2+2*(-1+x)^2*x*(x+2*(-1+x)*y0[1])*(z[3])^2))/x^2
        end
    else
        z[1]=0#?
    end
    
    
    return z[:]
end


# Defining the function in the RHS of the ution equation system
function SF_RHS(data,t,X)
    
    L=length(X)
    dy=zeros((L,length(data[1,:])));


    derxi_func = Spline1D(X[4:L-3], data[4:L-3,4],  k=4);

    funcs=[derxi_func];



    # update m, beta and psi data
    y0=[0.0 0.0 0.0]
    if twod==true
        data[4:L-3,1:3] = twod_n_rk4wrapper(RHS,y0,X[4:L-3],t,funcs,data[:,:])
    else
        data[4:L-3,1:3] = n_rk4wrapper(RHS,y0,X[4:L-3],t,funcs,data[:,:])
    end
    data=ghost(data)

    if twod==true
        epsilon=0.02
    else
        epsilon=0.0065
    end

    
        
    Threads.@threads for i in 4:L-3 #ORI
        
        dy[i,4]=derxi_evol(data,i,X) - dissipation(data,i,epsilon,4)[4]
        
    end
    dy[4,4]=extrapolate_in(dy[5,4], dy[6,4], dy[7,4], dy[8,4])

    dy[L-3,4]=extrapolate_out(dy[L-7,4], dy[L-6,4], dy[L-5,4], dy[L-4,4])

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

function masslossfunc(y)
    z=zeros(L)
    T00w=zeros(L)
    T01w=zeros(L)
    T11w=zeros(L)

    for i in 4:L-3

        if compactified==true
            r=y[i,8]
            der=y[i,5] #(x\[Chi]^(0,1))[u,r]=y[i,5]
        else
            x=y[i,8]
            r=x/(1-x)
            der=y[i,5]*(x-1)^2 #(x\[Chi]^(0,1))[u,r]=y[i,5]*(x-1)^2
        end
        
        T00w[i]=(1/(2*r^5*(1+r^2)^4))*(4*r^3*(1+r)^4*(1+r^2)^2*y[i,4]^2+exp(4*y[i,2])*(r-2*y[i,1])*((1+r^2)^2-(1+r^2+r*(1+r)*y[i,7])^2)^2+2*exp(2*y[i,2])*r^2*(1+r)^2*(1+r^2)*(r-2*y[i,1])*y[i,4]*((-1-2*r+r^2)*y[i,7]-r*(1+r+r^2+r^3)*der)+2*exp(4*y[i,2])*r*(r-2*y[i,1])^2*((-1-2*r+r^2)*y[i,7]-r*(1+r+r^2+r^3)*der)^2)

        T01w[i]=(1/(2*(r+r^3)^4))*exp(2*y[i,2])*((1+r^2)^4-2*(1+r^2)^2*(1+r^2+r*(1+r)*y[i,7])^2+(1+r^2+r*(1+r)*y[i,7])^4+2*exp(-2*y[i,2])*r^2*(1+r)^2*(1+r^2)*y[i,4]*((-1+(-2+r)*r)*y[i,7]-r*(1+r)*(1+r^2)*der)+2*r*(r-2*y[i,1])*((1-(-2+r)*r)*y[i,7]+r*(1+r)*(1+r^2)*der)^2)
    
        T11w[i]=(2*((1+2*r-r^2)*y[i,7]+r*(1+r+r^2+r^3)*der)^2)/(r^2*(1+r^2)^4)

        #T00s[i]=(\[Psi]^(1,0))[u,r]^2/r^2+(E^(4 y[i,2]) (r-2 y[i,1]) (y[i,3]-r y[i,6]) ((r-2 y[i,1]) (y[i,3]-r y[i,6])+2 E^(-2 y[i,2]) r^2 (\[Psi]^(1,0))[u,r]))/(2 r^6)

        #T01s[i]=(E^(2 y[i,2]) (r-2 y[i,1]) (y[i,3]-r y[i,6])^2)/(2 r^5)

        #T11s[i]=(y[i,3]-r y[i,6])^2/r^4

        z[i]=-4*exp(-2*y[i,2])*pi*r^2*T00w[i]+8*pi*r*T01w[i]*(r-2*y[i,1])-2*exp(2*y[i,2])*pi*T11w[i]*(r-2*y[i,1])^2
    
    end
    
    
    return z
end
function timeevolution(state_array,finaltime,run)#(state_array,finaltime,dir,run)

    t=0.0
    T_array = [0.0]
    iter = 0
    k=0
    massloss=zeros(L)
    lastprint_time=0.0
    while t<finaltime#@TRACK

        iter = iter + 1
        
        #update time increment

        if criticality!=true
            global dt = update_dt(initX,state_array[:,1],state_array[:,2],dt,ginit)      
        end
        t = t + dt
        
        if iter%500==0
            println("\n\niteration ", iter, " dt is ", dt, ", t=", t, " speed is ", speed(initX, state_array[:,1], state_array[:,2]), ", dx/dt=", dx/dt)
        end
        #println("\n\niteration ", iter, " dt is ", dt, ", t=", t, " speed is ", speed(initX, state_array[:,1], state_array[:,2]), ", dx/dt=", dx/dt)
        
        
        T_array = vcat(T_array,t)

        X=initX #state_array[:,5]
        X1=X[4:L-3]

        #evolve psi,r and xchi
        if twod==true
            state_array[:,:] = twod_rungekutta4molstep(SF_RHS,state_array[:,:],T_array,iter,X) #evolve psi,x
        else
            state_array[:,:] = rungekutta4molstep(SF_RHS,state_array[:,:],T_array,iter,X) #evolve psi,x
        end
        
        
        derxi_func = Spline1D(initX[4:L-3], state_array[4:L-3,4],  k=4);

        funcs=[derxi_func];

        #evolve m and beta together, new
        y0=[0.0 0.0 0.0]
        
        if twod==true
            state_array[4:L-3,1:3] = twod_n_rk4wrapper(RHS,y0,X1,t,funcs,state_array[:,:])
        else
            state_array[4:L-3,1:3] = n_rk4wrapper(RHS,y0,X1,t,funcs,state_array[:,:])
        end

        derderchi=Der_arrayLOP(state_array,4,X) .* (initX .- 1) .^ 2
        
        #threshold for apparent black hole formation
        if compactified==false
            global monitor_ratio[5:L-4] = 2 .* state_array[5:L-4,1] ./ initX[5:L-4]
        else
            global monitor_ratio[5:L-4] = 2 .* state_array[5:L-4,1] ./ initX[5:L-4] .* (1 .- initX[5:L-4])
        end

        run=int(run)
        if (iter%50==0)||((t>1.0)&&(t-lastprint_time)>0.01*(1.06-t))
            lastprint_time=t
            if zeroformat==true
                zero_print_muninn(files, t, [state_array[:,:] derderchi],res,"a")
            else
                print_muninn(files, t, [state_array[:,1:4] derderchi],res,"a",initX)
            end

        end

        if iter%50==0
            print_monitorratio("monitorratio", t, monitor_ratio[5:L-4],"a", initX[5:L-4])
        end
        """if hessian_control(state_array,t)==true
            global criticality = true
            global time = t
            break
        end"""

        

        
        if maximum(monitor_ratio)>0.71&&k==0
            global criticality = true
            k=k+1
            println("Supercritical evolution! At time ", t, ", iteration = ", iter)
            println("t = ", t, "iteration ", iter, " monitor ratio = ", maximum(monitor_ratio))
            global time = t
            break
        end

        
        if criticality == true
            break
        end
        
        if isnan(state_array[L-3,4])
            if criticality==false
                global time = t
                global explode = true
            end

            println("boom at time=", t)
            criticality=true
            break

        end

        #global time = t

        #println("state array at the end of iteration of time ", time, " is ", state_array)
        
    end
    
    if criticality==false
        global time = t
    end

    if t>1.4
        global time = 1.5
        global criticality = false
    end
    
    global evol_stats = [criticality A sigma r0 time explode run]

    return evol_stats, T_array

end