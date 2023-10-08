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

function init_gaussian_der(r,r0,sigma,A)
    n=length(r);
    if compactified==false
        if n==1
            z= A * (2 * exp(-(r-r0)^2/sigma^2) * r - 2 * exp(-(r-r0)^2/sigma^2) * (r-r0)*r^2/sigma^2)#exp(-((x/(1-x)-r0)/sigma)^2) * (3 * x^2 / (1-x) ^4 - (x/(1-x))^3 * (2*(x-r0*(-x+1)))/(sigma^2*(1-x)^3))
        else
            z=zeros(n);
            for i in 1:n
                rr = r[i]
                z[i] = A * (2 * exp(-(rr-r0)^2/sigma^2) * rr - 2 * exp(-(rr-r0)^2/sigma^2) * (rr-r0)*rr^2/sigma^2)
            end
        end
    else # inputted argument r is actually an x

        if n==1
            x=r
            r=x/(1-x)
            z= A * (2 * exp(-(r-r0)^2/sigma^2) * r - 2 * exp(-(r-r0)^2/sigma^2) * (r-r0)*r^2/sigma^2)#exp(-((x/(1-x)-r0)/sigma)^2) * (3 * x^2 / (1-x) ^4 - (x/(1-x))^3 * (2*(x-r0*(-x+1)))/(sigma^2*(1-x)^3))
        else
            z=zeros(n);
            for i in 1:n

                ## psi,r (corr?)
                x=r[i]
                rr = x/(1-x)
                z[i] = A * (2 * exp(-(rr-r0)^2/sigma^2) * rr - 2 * exp(-(rr-r0)^2/sigma^2) * (rr-r0)*rr^2/sigma^2)


                ## psi,x (correct)
                #x=r[i]
                #rr = x/(1-x)
                #z[i] = A * (2 * exp(-(rr-r0)^2/sigma^2) * rr - 2 * exp(-(rr-r0)^2/sigma^2) * (rr-r0)*rr^2/sigma^2) / (1-x)^2
                
            end
            z[n] = 0
        end
    end

    return z
end
"""
function new_init_gaussian_der(x,r0,sigma,A)
    n=length(x);
    if n==1
        z = A * exp(-((x/(1-x)-r0)/sigma)^2) * x^2/((x-1)^6*sigma^2) * (3*sigma^2 - 6*x*sigma^2 + x^2*(-2+3*sigma^2) - 2*r0*(x-1)*x)
        #z= A * exp(-((x/(1-x)-r0)/sigma)^2) * (3 * x^2 / (1-x) ^4 - (x/(1-x))^3 * (2*(x-r0*(-x+1)))/(sigma^2*(1-x)^3))
    else
        z=zeros(n);
        for i in 1:n
            if i<n-3 #avoid NaN for x=1, otherwise, it's 0
                #z[i] = 2*A * x[i]/(1-x[i])^3 * exp(-((x[i]/(1-x[i])-r0)/sigma)^2) * (1 - x[i]/(1-x[i]) * (x[i]/(1-x[i]) - r0) / sigma^2)
                z[i]=A * exp(-((x[i]/(1-x[i])-r0)/sigma)^2) * (3 * x[i]^2 / (1-x[i]) ^4 - (x[i]/(1-x[i]))^3 * (2*(x[i]-r0*(-x[i]+1)))/(sigma^2*(1-x[i])^3))
            end
        end
    end
    return z
end"""

function create_range(ori,stop,dx,N)

    array = []
    array = ori

    x=ori
    for i in 1:N
        x = x+dx
        array = vcat(array, x)
        
    end
    return array

end

# outputs xtilde(x)
function gridfunc(x)

    return 1/2 .+ 1/2 .* cos.( pi .* (1 .- 0.9 .* x)) #option 6
    
end;

# outputs dxtilde/dx(x)
function analytic_jacobian(x)

    if length(x) == 1
        if abs.(x .- 1.0)<10^(-15) #right
            println("hallo at analytic jacobian")
            return 0.0
        else
            return 0.9/2 * sin(pi * (1 - 0.9 * x)) #option6
        end
    else
        z = zeros(length(x))
        for i in 1:length(x)
            
            z[i] = 0.9 ./ 2 .* sin.(pi .* (1 .- 0.9 .* x[i])) #option6
        end
        
        return z
    end
end;

function inverse(x)
    return x = 1/2+1/2*cos(pi*(1-x))#1-1/pi*acos(2x-1)#-(acos(2*x-1)-pi)/(0.9*pi)
end

# Interpolation

function interpolate(x,x1,x2,y1,y2)
    return y1+(y2-y1)*(x-x1)/(x2-x1)
end


# Extrapolation

function extrapolate_out(y0,y1,y2,y3)
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
    g=abs.((1.0 .- initX[ori+1:L-4]) .^ 3.0 .* exp.(2 .* state_array[ori+1:L-4,2]) .* (2 .* state_array[ori+1:L-4,1] .- initX[ori+1:L-4] ./ (1 .- initX[ori+1:L-4])) ./ (2 .* initX[ori+1:L-4]))
    #println(g)
    
    #println(g)
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

    """if dt*(ginit/g) < 0.00014
        println("dt ", dt "\n g", g, " ginit ", ginit)
    end"""

    return  dx/g*0.5#dt*(ginit/g)

end

function find_origin(X)

    origin_i=Int64
    L=length(X)
    #origin_i = 4
    for i in 1:L
        if X[i] >= 0
            origin_i = i
            break
        end
    end
    #println("\nOrigin of the grid is at t = ", t, " is i = ", origin_i, "\n")

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

function rungekutta4(f,y0,T,func)
    n = length(T)
    y = zeros(n)
    y[1] = y0;
    for i in 1:n-1
        h = T[i+1] .- T[i]
        k1 = f(y[i], T[i],func)
        k2 = f(y[i] .+ k1 * h/2, T[i] .+ h/2,func)
        k3 = f(y[i] .+ k2 * h/2, T[i] .+ h/2,func)
        k4 = f(y[i] .+ k3 * h, T[i] .+ h,func)
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

function integrator(X,derpsi_func,data)

    L=length(X)
    ori=find_origin(X)

    integral = zeros(L)

    #taylor
    auxdata=zeros(L,4)
    auxdata[4:L-3,4]=DDer_array(state_array,4,initX)
    
    D3phi = auxdata[4,4]

    auxdata2=zeros(L,4)
    auxdata2[4:L-3,4]=DDer_array(auxdata,4,initX)
    D5phi = auxdata2[4,4]
    


    for i in ori:L-3

        if i == 4
            integral[i] = 0.0
        elseif i == 5
            h = X[i]-X[i-1]
            #integral[i] = h * (55/24*derpsi_func(X[i])-59/24*derpsi_func(X[i+1])+37/24*derpsi_func(X[i+2])-9/24*derpsi_func(X[i+3]))
            integral[i] = data[4,4] + 3*D3phi*X[i]^2/(3*2) + 5*D5phi*X[i]^4/(5*4*3*2)
        elseif i == 6
            h = X[i] - X[i-2]
            #integral[i] = h * (55/24*derpsi_func(X[i])-59/24*derpsi_func(X[i+1])+37/24*derpsi_func(X[i+2])-9/24*derpsi_func(X[i+3]))
            integral[i] = data[4,4] + 3*D3phi*X[i]^2/(3*2) + 5*D5phi*X[i]^4/(5*4*3*2)
        elseif i == 7
            h = X[i] - X[i-3]
            integral[i] = data[4,4] + 3*D3phi*X[i]^2/(3*2) + 5*D5phi*X[i]^4/(5*4*3*2)
            #integral[i] = h * (55/24*derpsi_func(X[i])-59/24*derpsi_func(X[i+1])+37/24*derpsi_func(X[i+2])-9/24*derpsi_func(X[i+3]))
        else
            h = X[i]-X[i-1]
            integral[i] = integral[i-1] + h * (55/24*derpsi_func(X[i-1]) - 59/24*derpsi_func(X[i-2]) + 37/24*derpsi_func(X[i-3]) - 9/24*derpsi_func(X[i-4]))
        end
    end

    return integral
end


function print_muninn(files, t, data, res, mode)
    #mode is "a" for append or "w" for write
    j=1
    if bisection==false
        for fl in files #normal run
            
            open(dir*"/muninnDATA/res$res/$fl.txt", mode) do file
                @printf file "\"Time = %.10e\n" t
                for i in 1:length(data[:,1])
                    @printf file "% .10e % .10e\n" data[i,5] data[i,j]
                end
                println(file) # insert empty line to indicate end of data set
                end
            j=j+1
        end
    else
        for fl in files #bisection search
            
            open(dir*"/muninnDATA/run$run/$fl.txt", mode) do file
                @printf file "\"Time = %.10e\n" t
                for i in 1:length(data[:,1])
                    @printf file "% .10e % .10e\n" data[i,5] data[i,j]
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
    y[3,1:4]=extrapolate_in(y[4,1:4], y[5,1:4], y[6,1:4], y[7,1:4])
    y[2,1:4]=extrapolate_in(y[3,1:4], y[4,1:4], y[5,1:4], y[6,1:4])
    y[1,1:4]=extrapolate_in(y[2,1:4], y[3,1:4], y[4,1:4], y[5,1:4])

    #outer boundary extrapolation
    y[L-2,1:4]=extrapolate_out(y[L-6,1:4], y[L-5,1:4], y[L-4,1:4], y[L-3,1:4])
    y[L-1,1:4]=extrapolate_out(y[L-5,1:4], y[L-4,1:4], y[L-3,1:4], y[L-2,1:4])
    y[L,1:4]=extrapolate_out(y[L-4,1:4], y[L-3,1:4], y[L-2,1:4], y[L-1,1:4])
   
    return y
end
function dissipation6(y,i,eps)
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
function dissipation4(y,i,eps)#0.02
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

"""function dissipation4(y,i,eps)#0.02
    if i==4
        delta4=(-13/6*y[i+5,:]+71/6*y[i+4,:]-77/3*y[i+3,:]+83/3*y[i+2,:]-89/6*y[i+1,:]+19/6*y[i,:])
    elseif i==5
        delta4=(-7/6*y[i+4,:]+41/6*y[i+3,:]-47/3*y[i+2,:]+53/3*y[i+1,:]-59/6*y[i,:]+13/6*y[i-1,:])
    elseif i==6 || i==7
        delta4=(y[i+2,:]-4*y[i+1,:]+6*y[i,:]-4*y[i-1,:]+y[i-2,:]);
    elseif i==8
        delta4=(-7/6*y[i-4,:]+41/6*y[i-3,:]-47/3*y[i-2,:]+53/3*y[i-1,:]-59/6*y[i,:]+13/6*y[i+1,:])
    else
        delta4=(-13/6*y[i-5,:]+71/6*y[i-4,:]-77/3*y[i-3,:]+83/3*y[i-2,:]-89/6*y[i-1,:]+19/6*y[i,:])        
    end

    return (-1)^2*eps*1/(dx)*delta4
end"""


#2nd order  dissipation, added to 1st order scheme
function dissipation2(y,i)
        delta2=(y[i+1,:]-2*y[i,:]+y[i-1,:]);
    return (-1)^1*epsilon*1/(dx)*delta2
end

"""#calculate dx~/dx at known gridpoints
function jacobian(y,X,i)
    if i==4
        z = (y[i+3]-4*y[i+2]+7*y[i+1]-4*y[i])/(2*(X[i+1]-X[i]))
    elseif i==L-3
        z = (-y[i-3]+4*y[i-2]-7*y[i-1]+4*y[i])/(2*(X[i]-X[i-1]))
    else
        z = (y[i+1]-y[i-1])/(2*(X[i+1]-X[i]))
    end

    return z
end

function jacobian(y,X)
    z=zeros(length(y))

    for i in 4:L-3
        if i==4
            z[i] = (y[i+3]-4*y[i+2]+7*y[i+1]-4*y[i])/(2*(X[i+1]-X[i]))
        elseif i==L-3
            z[i] = (-y[i-3]+4*y[i-2]-7*y[i-1]+4*y[i])/(2*(X[i]-X[i-1]))
        elseif i<L-3||i>4
            z[i] = (y[i+1]-y[i-1])/(2*(X[i+1]-X[i]))
        end
    end

    return z
end

#calculate dx~/dx at unknown gridpoints
function der_grid(X)
    
    aux = jacobian(X,originalX)
    
    dergrid_func = Spline1D(X[4:L-3],aux[4:L-3],k=4, bc="extrapolate")

    return dergrid_func
end"""

# Finite difference approximation
"""function Der(y,i,k,X)

    jacob = 1.0
    "if loggrid==true
        X = originalX
        jacob = jacobian_func(X[i])
    end"

    if i==4 # left boundary LOP1, TEM
        z = (y[i+3,k]-4*y[i+2,k]+7*y[i+1,k]-4*y[i,k])/(2*(X[i+1]-X[i]))*jacob
    elseif i==L-3
        z = (-y[i-3,k]+4*y[i-2,k]-7*y[i-1,k]+4*y[i,k])/(2*(X[i]-X[i-1]))*jacob
    else
        z = (y[i+1,k]-y[i-1,k])/(2*(X[i+1]-X[i]))*jacob
    end
        
    return z
    
end"""

function Der(y,i,k,X)

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
    
    return z
end

# Finite difference approximation
function Dertest(y,i,X)

    jacob = 1.0
    """if loggrid==true
        X = originalX
        jacob = jacobian_func(X[i])
    end"""
    
    if i==4 # left boundary LOP1, TEM
        z = (y[i+3]-4*y[i+2]+7*y[i+1]-4*y[i])/(2*(X[i+1]-X[i]))*jacob
    elseif i==L-3
        z = (-y[i-3]+4*y[i-2]-7*y[i-1]+4*y[i])/(2*(X[i]-X[i-1]))*jacob
    else
        z = (y[i+1]-y[i-1])/(2*(X[i+1]-X[i]))*jacob
    end
        
    return z
    
end


# Finite difference approximation
"""function DDer(y,i,k,X) #4th

    if i==4 # left boundary LOP2, TEM
        z = (15/4*y[i,k]-77/6*y[i+1,k]+107/6*y[i+2,k]-13*y[i+3,k]+61/12*y[i+4,k]-5/6*y[i+5,k])/((X[i+1]-X[i]))
        #z = (54*y[i+6,k]-208*y[i+5,k]+349*y[i+4,k]-336*y[i+3,k]+196*y[i+2,k]-64*y[i+1,k]+9*y[i,k])/(12*(X[i+1]-X[i])) #THIS NEEDS TO BE FIXED
    elseif i==5 # left boundary LOP1, TEM
        z = (-y[i+5,k]+7*y[i+4,k]-21*y[i+3,k]+34*y[i+2,k]-19*y[i+1,k]-9*y[i,k]+9*y[i-1,k])/(12*(X[i+1]-X[i]))
    else # central
        z = (-y[i+2,k]+16*y[i+1,k]-30*y[i,k]+16*y[i-1,k]-y[i-2,k])/(12*(X[i+1]-X[i]))
    

    end
    
    return z
end"""


#matrix
function unevenDer(y,i,k,X,spls)

    if k==4 #array of spline has variables m, beta and derpsi
        spl=spls[3]
    else
        spl=spls[k]
    end

    #for L-3 AND L-4 USE CENTRAL 
    

    dx=X[i+1]-X[i] #shouldnt this dx be constant, for error of derivatives to match?

    if i==4 # left boundary LOP1, TEM
        dx=X[i+1]-X[i]
        z = (-27*y[i,k]+58*y[i+1,k]-56*spl(X[i]+2*dx)+36*spl(X[i]+3*dx)-13*spl(X[i]+4*dx)+2*spl(X[i]+5*dx))/(12*(dx))
    elseif i==5 # left boundary LOP1, TEM
        dx=X[i+1]-X[i]
        z = (-2*spl(X[i]-dx)-15*y[i,k]+28*y[i+1,k]-16*spl(X[i]+2*dx)+6*spl(X[i]+3*dx)-spl(X[i]+4*dx))/(12*(dx))
    """elseif i==L-3
        dx=X[i]-X[i-1]
        z = -(-27*y[i,k]+58*y[i-1,k]-56*spl(X[i]-2*dx)+36*spl(X[i]-3*dx)-13*spl(X[i]-4*dx)+2*spl(X[i]-5*dx))/(12*(dx)) #12-06 i did *-1
    elseif i==L-4 # right boundary TEM
        dx=X[i+1]-X[i]
        z = -(-2*y[i+1,k]-15*y[i]+28*spl(X[i]-dx)-16*spl(X[i]-2*dx)+6*spl(X[i]-3*dx)-spl(X[i]-4*dx))/(12*(dx)) #12-06 i did *-1"""
    else
        

        if(X[i]-dx)<0.0||(X[i]-2*dx)<0.0 #avoid evaluating spline out of domain
            dx=X[i+1]-X[i]
            z = (-27*y[i,k]+58*y[i+1,k]-56*spl(X[i]+2*dx)+36*spl(X[i]+3*dx)-13*spl(X[i]+4*dx)+2*spl(X[i]+5*dx))/(12*(dx))
            #println("warnign at ori, i is ", i)
        """elseif(X[i]+dx)>1.0||(X[i]+2*dx)>1.0 #avoid evaluating spline out of domain
            dx=X[i]-X[i-1]
            #println("warning!, X[i]", X[i], "i is ", i)
            #z = (spl(X[i]+3*dx)-4*spl(X[i]+2*dx)+7*y[i+1,k]-4*y[i,k])/(2*dx)
            z = -(-27*y[i,k]+58*y[i-1,k]-56*spl(X[i]-2*dx)+36*spl(X[i]-3*dx)-13*spl(X[i]-4*dx)+2*spl(X[i]-5*dx))/(12*(dx))"""
        else
            dx=X[i+1]-X[i]
            z = (-spl(X[i]+2dx)+8*y[i+1,k]-8*spl(X[i]-dx)+spl(X[i]-2*dx))/(12*(dx))
        end
    end
        
    return z
    
end




function DDer_array(y,k,X)

    aux=zeros(L-6)
    i=1
    j=4
    for xx in X[4:L-3]
        aux[i] = Der(y,j,k,X)
        i=i+1
        j=j+1
    end

    auxdata=zeros(L,4)
    auxdata[4:L-3,k]=aux[:]
    aux2=zeros(L-6)
    i=1
    j=4
    for xx in X[4:L-3]
        aux2[i] = Der(auxdata,j,k,X)
        i=i+1
        j=j+1
    end
    
    return aux2
end

int(x) = floor(Int, x)

function chebyshev(N)

    X=zeros(N)
    
    for i in 1:N
        #X[i]=1/2+1/2*cos((2*i-1)*pi/(2*N))
        if i==1
            X[i]=0.0
        else
            X[i]=(1/2+1/2*cos((2*i-1)*pi/(2*N)))
        end
    end

    return sort(X)
end

function chebyshev_weigth(X)
    w=ones(length(X))
    len=length(X)
    for i in 1:len
            
        w[i]=w[i]=1/2+1/2*abs(cos((i-1)*pi/(len)))#1/2+1/2*(cos(1/2*(i-1)*pi/(len)))

    end
    return w
end
function chebyshev_cut(X)
    N=length(X)
    new_grid=zeros(int(N/4))
    
    new_grid[1:int(N/4)] = X[1:int(N/4)]
    new_grid=vcat(new_grid, X[int(N/4):4:int(3*N/4)])
    new_grid=vcat(new_grid, X[int(3*N/4):2:int(N)])
    
    #deleteat!(A, 2)
    return new_grid
end

function h(z)
    return acos(-1 + 2*z)
end

function bulkSF(y,i,X)
    
    #psi,x
    if compactified == false
        r=X[i]
        dy=(1/(2*r^3))*exp(2*y[i,2])*(-2*y[i,1]*y[i,3]+2*r*y[i,3]*Der(y,i,1,X)-2*r^2*y[i,3]*Der(y,i,2,X)+4*r*y[i,1]*y[i,3]*Der(y,i,2,X)+2*r*y[i,1]*y[i,4]-2*r^2*Der(y,i,1,X)*y[i,4]+2*r^3*Der(y,i,2,X)*y[i,4]-4*r^2*y[i,1]*Der(y,i,2,X)*y[i,4]+r^3*Der(y,i,4,X)-2*r^2*y[i,1]*Der(y,i,4,X))
    else
        ## psi,r evol equation
        x=X[i]
        r = x/(1-x)
        dy=(1/(2*r^3))*exp(2*y[i,2])*(r^3*(1-x)^2*Der(y,i,4,X)-2*r^2*(1-x)^2*Der(y,i,4,X)*y[i,1]+2*r*y[i,1]*y[i,4]-2*y[i,1]*y[i,3]-2*r^2*(1-x)^2*y[i,4]*Der(y,i,1,X)+2*r*(1-x)^2*y[i,3]*Der(y,i,1,X)+2*r^3*(1-x)^2*y[i,4]*Der(y,i,2,X)-4*r^2*(1-x)^2*y[i,1]*y[i,4]*Der(y,i,2,X)-2*r^2*(1-x)^2*y[i,3]*Der(y,i,2,X)+4*r*(1-x)^2*y[i,1]*y[i,3]*Der(y,i,2,X))



        ## psi,x evolu equation
        #x=X[i]
        #dy=(1/(2*x^3))*exp(2*y[i,2])*(2*(-1+x)*y[i,1]*((-1+x)^2*y[i,3]*(1+2*(-1+x)*x*Der(y,i,2,X))+x*((-1+x)^3*(1+2*(-1+x)*x*Der(y,i,2,X))*y[i,4]+x*Der(y,i,4,X)))+x*(2*(-1+x)^3*y[i,3]*((-1+x)*Der(y,i,1,X)+x*Der(y,i,2,X))+x*(2*(-1+x)^5*Der(y,i,1,X)*y[i,4]+x*(2*(-1+x)^4*Der(y,i,2,X)*y[i,4]+Der(y,i,4,X)))))
        
        """if loggrid==false
            x=X[i]
            dy=(1/(2*x^3))*exp(2*y[i,2])*(2*(-1+x)*y[i,1]*((-1+x)^2*y[i,3]*(1+2*(-1+x)*x*Der(y,i,2,X))+x*((-1+x)^3*(1+2*(-1+x)*x*Der(y,i,2,X))*y[i,4]+x*Der(y,i,4,X)))+x*(2*(-1+x)^3*y[i,3]*((-1+x)*Der(y,i,1,X)+x*Der(y,i,2,X))+x*(2*(-1+x)^5*Der(y,i,1,X)*y[i,4]+x*(2*(-1+x)^4*Der(y,i,2,X)*y[i,4]+Der(y,i,4,X)))))
        else
            dy=(1/(4*pi^2*(pi-h(x))^3))*exp(2*y[i,2])*(-4*h(x)*y[i,1]*(pi*h(x)^2*y[i,3]*(pi-2*sqrt(-((-1+x)*x))*(pi-h(x))*h(x)*Der(y,i,2,X))-(pi-h(x))*(pi*sqrt(-((-1+x)*x))*h(x)^3*y[i,4]+2*(-1+x)*x*(pi-h(x))*h(x)^4*Der(y,i,2,X)*y[i,4]-pi^2*(pi-h(x))*Der(y,i,4,X)))+(pi-h(x))*(4*pi*sqrt(-((-1+x)*x))*h(x)^3*y[i,3]*(h(x)*Der(y,i,1,X)-(pi-h(x))*Der(y,i,2,X))-(pi-h(x))*(-4*(-1+xt)*x*h(x)^5*Der(y,i,1,X)*y[i,4]-(pi-h(x))*(-4*(-1+xt)*x*h(x)^4*Der(y,i,2,X)*y[i,4]+2*pi^2*Der(y,i,4,X)))))

        end"""
    end
    return dy
end


function bulkSF(y,i,X,spls)
    
    #psi,x
    if compactified == false
        r=X[i]
        dy=(1/(2*r^3))*exp(2*y[i,2])*(-2*y[i,1]*y[i,3]+2*r*y[i,3]*Der(y,i,1,X)-2*r^2*y[i,3]*Der(y,i,2,X)+4*r*y[i,1]*y[i,3]*Der(y,i,2,X)+2*r*y[i,1]*y[i,4]-2*r^2*Der(y,i,1,X)*y[i,4]+2*r^3*Der(y,i,2,X)*y[i,4]-4*r^2*y[i,1]*Der(y,i,2,X)*y[i,4]+r^3*Der(y,i,4,X)-2*r^2*y[i,1]*Der(y,i,4,X))
    else
        ## psi,r evol equation
        x=X[i]
        r = x/(1-x)
        dy=(1/(2*r^3))*exp(2*y[i,2])*(r^3*(1-x)^2*unevenDer(y,i,4,X,spls)-2*r^2*(1-x)^2*unevenDer(y,i,4,X,spls)*y[i,1]+2*r*y[i,1]*y[i,4]-2*y[i,1]*y[i,3]-2*r^2*(1-x)^2*y[i,4]*unevenDer(y,i,1,X,spls)+2*r*(1-x)^2*y[i,3]*unevenDer(y,i,1,X,spls)+2*r^3*(1-x)^2*y[i,4]*unevenDer(y,i,2,X,spls)-4*r^2*(1-x)^2*y[i,1]*y[i,4]*unevenDer(y,i,2,X,spls)-2*r^2*(1-x)^2*y[i,3]*unevenDer(y,i,2,X,spls)+4*r*(1-x)^2*y[i,1]*y[i,3]*unevenDer(y,i,2,X,spls))


    end
    return dy
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

function h(x)
    return acos(-1+2*x)
end

#EXACTLY THE SAME POINTS MUST BE CALCULATED THE SAME WAY
function RHS(y0,x1,time,func,i,data)
    
    z=zeros(length(y0))
    derpsi = func

    if compactified==false
        r=x1
        z[3]=derpsi(r)
    else
        #psi,x
        z[3]=derpsi(x1)/(1-x1)^2
    end
    
    
    
    #m and beta
    if x1<10^(-15) #left
        z[1] = 0.0;
        z[2] = 0.0;
    else
    
        if compactified == false
            r=x1
            z[1] = (r - 2.0 * y0[1]) * 2.0 .* pi .* r * ((r*z[3]-y0[3])/r^2.0) ^ 2.0
            z[2] = 2.0 .* pi .* r * ((r*z[3]-y0[3])/r^2.0) ^ 2.0
        else
            if loggrid==false
                x=x1
                z[1] = - 2.0 .* pi .* (-1.0 .+ x) .* (y0[3] .+ (-1 + x) .* x .* z[3]) .^ 2.0 ./ x .^ 3.0 .* ( x ./ (1.0 .-x ) .- 2 .* y0[1])
                z[2] = - 2.0 .* pi .* (-1.0 .+ x) .* (y0[3] .+ (-1 + x) .* x .* z[3]) .^ 2.0 ./ x .^ 3.0
                if abs.(x1 .- 1.0)<10^(-15)
                    z[1] = 0.0
                    z[2] = 0.0
                    z[3] = 0.0
                end
            else
                #x = x1
                #z[1] = (2.0 .* h(x) .* (pi .* y0[3] + sqrt(-((-1+x) .* x)) .* h(x) .* (-pi .+ h(x)) .* z[3])^2)/(sqrt(-((-1+x) .* x)) .* (pi-h(x))^3) .* ((pi - h(x) .* (1.0 .+ 2.0 .* y0[1])))/h(x)
                #z[2] = (2.0 .* h(x) .* (pi .* y0[3] + sqrt(-((-1+x) .* x)) .* h(x) .* (-pi .+ h(x)) .* z[3])^2)/(sqrt(-((-1+x) .* x)) .* (pi-h(x))^3)
                x=x1
                z[1] = - 2.0 .* pi .* (-1.0 .+ x) .* (y0[3] .+ (-1 + x) .* x .* z[3]) .^ 2.0 ./ x .^ 3.0 .* ( x ./ (1.0 .-x ) .- 2 .* y0[1])
                z[2] = - 2.0 .* pi .* (-1.0 .+ x) .* (y0[3] .+ (-1 + x) .* x .* z[3]) .^ 2.0 ./ x .^ 3.0
                if abs.(x1 .- 1.0)<10^(-15)
                    z[1] = 0.0
                    z[2] = 0.0
                    z[3] = 0.0
                end
            end
        end

    end

    return z[:]
end


# Defining the function in the RHS of the evolution equation system
using Base.Threads

function SF_RHS(data,t,X)
    
    L=length(X)
    dy=zeros((L,length(data[1,:])));

    # update interpolation of psi,x
    derpsi_func = Spline1D(X[4:L-3],data[4:L-3,4],k=4)
    
    # update m, beta and psi data
    y0=[0.0 0.0 0.0]
    data[4:L-3,1:3] = n_rk4wrapper(RHS,y0,X[4:L-3],t,derpsi_func,data[:,:])
    
    #NEW
    if loggrid==true
        m_func = Spline1D(X[4:L],data[4:L,1],k=4)
        beta_func = Spline1D(X[4:L],data[4:L,2],k=4)

        funcs=[m_func beta_func derpsi_func]
    end

    if loggrid==false
        Threads.@threads for i in 4:L-3 #ORI
            if X[i]<10^(-15) #left
                dy[i,4]= +1.0/2.0 * (1/(1-X[i])^2 * Der(data,i,4,X) + 2/(1-X[i])^3*data[i,4])  - dissipation6(data,i,0.0065)[4]#- dissipation4(data,i,0.02)[4];
                
            elseif abs.(X[i] .- 1.0)<10^(-15)
                dy[i,4]= 0.0 - dissipation6(data,i,0.0065)[4]#- dissipation4(data,i,0.02)[4]
                
            else
                dy[i,4]=bulkSF(data,i,X) - dissipation6(data,i,0.0065)[4]#- dissipation4(data,i,0.02)[4]
            end
        end
    else
        Threads.@threads for i in 4:L-3 #ORI
            if X[i]<10^(-15) #left
                dy[i,4]= +1.0/2.0 * (1/(1-X[i])^2 * unevenDer(data,i,4,X,funcs) + 2/(1-X[i])^3*data[i,4])  - dissipation6(data,i,0.0065)[4]#- dissipation4(data,i,0.02)[4];
                
            elseif abs.(X[i] .- 1.0)<10^(-15)
                dy[i,4]= 0.0 - dissipation6(data,i,0.0065)[4]#- dissipation4(data,i,0.02)[4]
                
            else
                dy[i,4]=bulkSF(data,i,X,funcs) - dissipation6(data,i,0.0065)[4]#- dissipation4(data,i,0.02)[4]
            end
        end
    end
    
    
    dy[4,4]=extrapolate_in(dy[5,4], dy[6,4], dy[7,4], dy[8,4])
  
    #dy=ghost(dy)
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
    
    dy=ghost(dy)
    return dy
end

#double the grid and get the data there
function doublegrid(y)
    
    L1=length(y[:,5])
    newX = y[1:3,5]
    newX = vcat(newX,doubleX(y[4:L1-3,5]))
    newX =vcat(newX,y[L1-2:L1,5])
    
    m_func = Spline1D(y[4:L1-3,5],y[4:L1-3,1],k=4)
    beta_func = Spline1D(y[4:L1-3,5],y[4:L1-3,2],k=4)
    psi_func = Spline1D(y[4:L1-3,5],y[4:L1-3,3],k=4)
    derspi_func = Spline1D(y[4:L1-3,5],y[4:L1-3,4],k=4)

    L = length(newX)

    z = zeros(L,5)

    z[:,5] = newX
    z[:,1] = m_func(newX)
    z[:,2] = beta_func(newX)
    z[:,3] = psi_func(newX)
    z[:,4] = derpsi_func(newX)

    return z

end

# double only the grid
function doubleX(X)
    new_grid=[X[1]]

    L = length(X)

    for i in 1:(L-1)
        h = X[i+1]-X[i]
        
        new_grid = vcat(new_grid, X[i]+h/2)
        new_grid = vcat(new_grid, X[i+1])
    end

    return new_grid
end

#using ProgressMeter

function timeevolution(state_array,finaltime,dir,run)

    t=0.0
    T_array = [0.0]
    iter = 0
    

    while t<finaltime#@TRACK

        iter = iter + 1

        #update time increment
        #global dt = update_dt(initX,state_array[:,1],state_array[:,2],dt,ginit)
        #global dt=0.0000000001
        t = t + dt
        if iter%20==0
            println("\n\niteration ", iter, " dt is ", dt, ", t=", t, " speed is ", speed(initX, state_array[:,1], state_array[:,2]), ", dx/dt=", dx/dt)
        end
        #println("\n\niteration ", iter, " dt is ", dt, ", t=", t, " speed is ", speed(initX, state_array[:,1], state_array[:,2]), ", dx/dt=", dx/dt)

        T_array = vcat(T_array,t)

        #X = update_grid(state_array[:,:],T,t)
        
        X=initX #state_array[:,5]
        X1=X[4:L-3]
       
        #evolve psi,x
        state_array[:,:] = rungekutta4molstep(SF_RHS,state_array[:,:],T_array,iter,X) #evolve psi,x AQUII
        state_array=ghost(state_array)
    
        # update interpolation of psi,x
        derpsi_func = Spline1D(X[4:L-3],state_array[4:L-3,4],k=4)#new

        #evolve m and beta together, new
        y0=[0.0 0.0 0.0]
        state_array[4:L-3,1:3] = n_rk4wrapper(RHS,y0,X1,t,derpsi_func,state_array[:,:])
        state_array=ghost(state_array)
        

        run=int(run)
        
        if iter%100==0
            print_muninn(files, t, state_array[:,1:5],res,"a")
        end
        #print_muninn(files, t, state_array[:,1:5],res,"a")

        #threshold for apparent black hole formation
        if compactified==false
            global monitor_ratio[5:L-4] = 2 .* state_array[5:L-4,1] ./ initX[5:L-4]
        else
            global monitor_ratio[5:L-4] = 2 .* state_array[5:L-4,1] ./ initX[5:L-4] .* (1 .- initX[5:L-4])
        end


       
        if maximum(monitor_ratio)>0.70
            global criticality = true
            println("Supercritical evolution! At time ", t, ", iteration = ", iter)
            println("t = ", t, "iteration ", iter, " monitor ratio = ", maximum(monitor_ratio))
            global time = t
        end



        """if criticality == true
            break
        end"""
        
        if isnan(state_array[L-3,4])
            global explode = true
            println("boom at time=", t, " timestep = ", iter)
            global timestep = iter
            break
        end

        global time = t
        
    end
    
    

    global evol_stats = [criticality A sigma r0 time explode run]

    return evol_stats, T_array

end    

function epsilon(X,i,dt,dx)
    #minimum([dx/dt*(1/2)^(2*3), 10])
    #println("dissipation epsilon is ", (dx/dt*(1/2)^(2*3)))
    if i != L-3
        dx=X[i+1]-X[i]
    elseif i==L-3
        dx = X[i]-X[i-1]
    end
    return (dx/dt*(1/2)^(2*3+1))
end

function epsilon(dt,dx)
    #minimum([dx/dt*(1/2)^(2*3), 10])
    #println("dissipation epsilon is ", (dx/dt*(1/2)^(2*3)))
    
    return (dx/dt*(1/2)^(2*3+1))
end

function twod_epsilon(dt,dx)

    
    return (dx/dt*(1/2)^(2*2+1))
end
