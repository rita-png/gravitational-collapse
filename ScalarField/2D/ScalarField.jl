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
    return dx/g*0.5
   

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
function print_muninn(files, t, data, res, mode)
    #mode is "a" for append or "w" for write
    j=1
    if bisection==false
        for fl in files #normal run
            
            open(dir*"/muninnDATA/res$res/$fl.txt", mode) do file
            #open("./DATA/muninnDATA/res$res/$fl.txt", mode) do file
                @printf file "\"Time = %.10e\n" t
                for i in 1:length(data[:,1])
                    @printf file "% .10e % .10e\n" data[i,5] data[i,j]
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


# Finite difference approximation
function Der(y,i,k,X)

    jacob = 1.0
"""    if loggrid==true
        X = originalX
        jacob = jacobian_func(X[i])
    end"""

    if i==4 # left boundary LOP1, TEM
        z = (y[i+3,k]-4*y[i+2,k]+7*y[i+1,k]-4*y[i,k])/(2*(X[i+1]-X[i]))*jacob
    elseif i==L-3
        z = (-y[i-3,k]+4*y[i-2,k]-7*y[i-1,k]+4*y[i,k])/(2*(X[i]-X[i-1]))*jacob
    else
        z = (y[i+1,k]-y[i-1,k])/(2*(X[i+1]-X[i]))*jacob
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
                x = x1
                z[1] = (2.0 .* h(x) .* (pi .* y0[3] + sqrt(-((-1+x) .* x)) .* h(x) .* (-pi .+ h(x)) .* z[3])^2)/(sqrt(-((-1+x) .* x)) .* (pi-h(x))^3) .* ((pi - h(x) .* (1.0 .+ 2.0 .* y0[1])))/h(x)
                z[2] = (2.0 .* h(x) .* (pi .* y0[3] + sqrt(-((-1+x) .* x)) .* h(x) .* (-pi .+ h(x)) .* z[3])^2)/(sqrt(-((-1+x) .* x)) .* (pi-h(x))^3)
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
    data[4:L-3,1:3] = twod_n_rk4wrapper(RHS,y0,X[4:L-3],t,derpsi_func,data[:,:])
    

    Threads.@threads for i in 4:L-3 #ORI
        if X[i]<10^(-15) #left
            dy[i,4]= 0.0 - dissipation4(data,i,0.02)[4];
            
        elseif abs.(X[i] .- 1.0)<10^(-15)
            dy[i,4]= 0.0 - dissipation4(data,i,0.02)[4]
            
        else
            dy[i,4]=bulkSF(data,i,X) - dissipation4(data,i,0.02)[4]
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


function timeevolution(state_array,finaltime,dir,run)

    t=0.0
    T_array = [0.0]
    iter = 0
    k = 0

    while t<finaltime#@TRACK

        iter = iter + 1

        #update time increment

        if criticality!=true#||dt>0.00000001
            global dt = update_dt(initX,state_array[:,1],state_array[:,2],dt,ginit)      
        end
        
        t = t + dt
        if iter%200==0
            println("\n\niteration ", iter, " dt is ", dt, ", t=", t, " speed is ", speed(initX, state_array[:,1], state_array[:,2]), ", dx/dt=", dx/dt)
        end
        #println("\n\niteration ", iter, " dt is ", dt, ", t=", t, " speed is ", speed(initX, state_array[:,1], state_array[:,2]), ", dx/dt=", dx/dt)

        T_array = vcat(T_array,t)

        #X = update_grid(state_array[:,:],T,t)
        
        X=initX #state_array[:,5]
        X1=X[4:L-3]
       
        #evolve psi,x
        state_array[:,:] = twod_rungekutta4molstep(SF_RHS,state_array[:,:],T_array,iter,X) #evolve psi,x
        state_array=ghost(state_array)
    
        # update interpolation of psi,x
        derpsi_func = Spline1D(X[4:L-3],state_array[4:L-3,4],k=4)#new

        #evolve m and beta together, new
        y0=[0.0 0.0 0.0]
        state_array[4:L-3,1:3] = twod_n_rk4wrapper(RHS,y0,X1,t,derpsi_func,state_array[:,:])
        state_array=ghost(state_array)
        

        run=int(run)

        if (iter%50==0&&t>0.3)||(t>0.85&&iter%10==0)
        #if iter%100==0
            print_muninn(files, t, state_array[:,1:5],res,"a")
        end
        #print_muninn(files, t, state_array[:,1:5],res,"a")

        #threshold for apparent black hole formation
        if compactified==false
            global monitor_ratio[5:L-4] = 2 .* state_array[5:L-4,1] ./ initX[5:L-4]
        else
            global monitor_ratio[5:L-4] = 2 .* state_array[5:L-4,1] ./ initX[5:L-4] .* (1 .- initX[5:L-4])
        end


       
        if maximum(monitor_ratio)>0.88 && k==0
            global criticality = true
            k=k+1
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
    
    if t>2.9
        global time = 3.0
        global criticality = false
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
