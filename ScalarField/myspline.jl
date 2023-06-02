###########################################################
##                  Spline interpolation                 ##
## k=3 working                                           ##
## k=4 needs one last constraint                         ##
## boundary constraints should be checked for both cases ##
###########################################################

using LinearSolve

function cubic_f(x::Float128,c::Array{Float128})


    #z = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4
    z = c[1]*x^3 + c[2]*x^2 + c[3]*x + c[4]

    return z

end

function quartic_f(x,c)

    #z = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4
    z = c[1]*x^4 + c[2]x^3 + c[3]*x^2 + c[4]*x + c[5]

    return z

end

function findInterval(x,vec::Array{Float64})
    out = zeros(Int,length(x))
    vec = unique(vec)
    #sort!(vec)

    for j in 1:length(x)
        if x[j] < vec[1]
            out[1] = -1
        elseif x[j] > vec[end]
            out[end] = -1
        else
            out[j] = searchsortedfirst(vec,x[j])-1
            if out[j] == 0
                out[j] = 1
            end
        end
    end
    return out
end


function evalInterval(x,vec::Array{Float128},sol,k)
    out = zeros(Int,length(x))
    vec = unique(vec)
    #sort!(vec)
    #z=zeros(length(x))
    z=Array{Float128}(undef, length(x))
    z.=0-0

    for j in 1:length(x)
        if x[j] < vec[1]
            out[1] = -1
        elseif x[j] > vec[end]
            out[end] = -1
        else
            out[j] = searchsortedfirst(vec,x[j])-1
            if out[j] == 0
                out[j] = 1
            end
        end
        z[j]=cubic_f(x[j],sol[(k+1)*(out[j]-1)+1:(k+1)*out[j]])#cubic_f(x[j],sol[(k+1)*(loc-1)+1:(k+1)*loc])
    end
    return z
end




function splinethree(x::Array{Float128},y::Array{Float128},x1::Array{Float128})

    N=length(x)
    k=3

    #A=zeros((k+1)*(N-1),(k+1)*(N-1))
    A=Matrix{Float128}(undef, (k+1)*(N-1),(k+1)*(N-1))
    A.=0

    #b=zeros((k+1)*(N-1))
    b=Array{Float128}(undef, (k+1)*(N-1))
    b.=0
    
    
    j=1
    #Si(xi)=yi
    for i in 1:N-1
        aux=[x[i]^3, x[i]^2, x[i], 1]
        
        A[j,(k+1)*(i-1)+1:(k+1)*i]=aux
        b[j]=y[i]
        j=j+1
    end

    """for i in 1:N-1
        aux=zeros((k+1)*(N-1))
        aux[(k+1)*(i-1)+1:(k+1)*i]=[x[i]^3, x[i]^2, x[i], 1]
        
        A[j,:]=aux
        b[j]=y[i]
        j=j+1
    end"""
 
    
    #Si(xi+1)=yi+1
    for i in 1:N-1
        aux=[x[i+1]^3, x[i+1]^2, x[i+1], 1]
        
        A[j,(k+1)*(i-1)+1:(k+1)*i]=aux
        b[j]=y[i+1]
        j=j+1
    end

    #S'i(xi+1)=S'i+1(xi+1)
    for i in 1:N-2
        aux=zeros((k+1)*(N-1))
        aux[(k+1)*(i-1)+1:(k+1)*(i+1)]=[3*x[i+1]^2, 2*x[i+1], 1, 0, -3*x[i+1]^2, -2*x[i+1], -1, 0]      
        A[j,:]=aux
        b[j]=0
        j=j+1
    end

    #S''i(xi+1)=S''i+1(xi+1)
    for i in 1:N-2
        aux=zeros((k+1)*(N-1))
        aux[(k+1)*(i-1)+1:(k+1)*(i+1)]=[6*x[i+1], 2, 0, 0, -6*x[i+1], -2, 0,0]
        A[j,:]=aux
        b[j]=0
        j=j+1
    end
    
    #endpoints constraint
    aux=zeros((k+1)*(N-1))
    aux[1:(k+1)*2]=[6*x[1], 2, 0, 0, 0, 0, 0,0]
    A[j,:]=aux
    b[j]=0
    j=j+1
    aux=zeros((k+1)*(N-1))
    #aux[(k+1)*(N-1)-7:(k+1)*(N-1)]=[0, 0, 0, 0, 6*x[N], 2, 0,0]
    aux[(k+1)*(N-1)-(k+1)*2+1:(k+1)*(N-1)]=[0, 0, 0, 0, 6*x[N], 2, 0,0]
    A[j,:]=aux
    b[j]=0
    
    #println(A)
    #println(" ")
    #println(b)
    prob = LinearProblem(A, b)
    sol = solve(prob)

"""    # evaluate the spline
    L=length(x1)
    z = zeros(L)
    
    #locs=findInterval(x1,x)
    
    "for j in 1:L
        loc=0
        for i in 1:N-1
            if x1[j]>x[i]&&x1[j]<=x[i+1]
                #println("x is between ",x[i], " and ", x[i+1])
                loc=i
                break
            elseif x1[j]==0
                loc=1
            end
        end
        
        if loc==0
            println("Warning: Trying to evaluate spline out of its domain!")
        end
    #z[j] = cubic_f(x1[j],sol[(k+1)*(loc-1)+1:(k+1)*loc])
    end"
    

    #println(cubic_f(x1,sol[(k+1)*(loc-1)+1:(k+1)*loc]))

    #println(cubic_f(x1,sol[(k+1)*(locs-1)+1:(k+1)*locs]))
    z=evalInterval(x1,x,sol,k)
    
    return z"""

    return sol
end

function splinefour(x,y,x1)

    N=length(x)
    k=4

    A=zeros((k+1)*(N-1),(k+1)*(N-1))
    b=zeros((k+1)*(N-1))

    j=1
    #Si(xi)=yi
    for i in 1:N-1
        aux=[x[i]^4, x[i]^3, x[i]^2, x[i], 1]
        
        A[j,(k+1)*(i-1)+1:(k+1)*i]=aux
        b[j]=y[i]
        j=j+1
    end


    #Si(xi+1)=yi+1
    for i in 1:N-1
        aux=[x[i+1]^4, x[i+1]^3, x[i+1]^2, x[i+1], 1]
        
        A[j,(k+1)*(i-1)+1:(k+1)*i]=aux
        b[j]=y[i+1]
        j=j+1
    end

    
    #S'i(xi+1)=S'i+1(xi+1)
    for i in 1:N-2
        aux=zeros((k+1)*(N-1))
        aux[(k+1)*(i-1)+1:(k+1)*(i+1)]=[4*x[i+1]^3, 3*x[i+1]^2, 2*x[i+1], 1, 0, -4*x[i+1]^3, -3*x[i+1]^2, -2*x[i+1], -1, 0]      
        A[j,:]=aux
        b[j]=0
        j=j+1
    end


    #S''i(xi+1)=S''i+1(xi+1)
    for i in 1:N-2
        aux=zeros((k+1)*(N-1))
        aux[(k+1)*(i-1)+1:(k+1)*(i+1)]=[12*x[i+1]^2, 6*x[i+1], 2, 0, 0, -12*x[i+1]^2, -6*x[i+1], -2, 0,0]
        A[j,:]=aux
        b[j]=0
        j=j+1
    end
    
    #endpoints constraint
    aux=zeros((k+1)*(N-1))
    aux[1:(k+1)*2]=[12*x[1]^2, 6*x[1], 2, 0, 0, 0, 0, 0, 0, 0]
    A[j,:]=aux
    b[j]=0
    j=j+1
    aux=zeros((k+1)*(N-1))
    aux[(k+1)*(N-1)-(k+1)*2+1:(k+1)*(N-1)]=[0, 0, 0, 0, 0, 12*x[N]^2, 6*x[N], 2, 0,0]
    A[j,:]=aux
    b[j]=0
    j=j+1
    aux=zeros((k+1)*(N-1))
    aux[1:(k+1)*2]=[24*x[1], 6, 2, 0, 0, 0, 0, 0, 0, 0]
    A[j,:]=aux
    b[j]=0
    j=j+1
    aux=zeros((k+1)*(N-1))
    aux[1:(k+1)*2]=[0 ,0, 0, 0, 0, 24*x[L], 6, 2, 0, 0]
    A[j,:]=aux
    b[j]=0

    println(A)

    #Missing one constraint for the 4th order spline
    """
    println(A[1:j,:])
    println(" ")
    println(b[1:j])
    prob = LinearProblem(A[1:j,:], b[1:j])
    sol = solve(prob)
    println(sol)
    
    loc=0
    for i in 1:N-1
        if x1>x[i]&&x1<x[i+1]
            println("x is between ",x[i], " and ", x[i+1])
            loc=i
        end
    end
    if loc==0
        println("Warning: Trying to evaluate spline out of its domain!")
    end

    println(quartic_f(x1,sol[(k+1)*(loc-1)+1:(k+1)*loc]))

    return sol[(k+1)*(loc-1)+1:(k+1)*loc]"""
end

"""

x=[0 1 2 4]
y=[0 2 4 8]
x1=[0.5 1.0 1.5]



solution=splinethree(x,y,x1)
println(solution)
#splinefour(x,y,1.5)


"""