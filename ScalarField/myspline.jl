

function cubic_f(x,c)

    #z = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4
    z = c[1]*x^3 + c[2]*x^2 + c[3]*x + c[4]

    return z

end

function quartic_f(x,c)

    #z = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4
    z = c[1]*x^4 + c[2]x^3 + c[3]*x^2 + c[4]*x + c[5]

    return z

end

function myspline(x,y,x1)

    N=length(x)

    A=zeros(4*(N-1),4*(N-1))
    b=zeros(4*(N-1))
    k=3

    j=1
    #Si(xi)=yi
    for i in 1:N-1
        aux=[x[i]^3, x[i]^2, x[i], 1]
        
        A[j,(k+1)*(i-1)+1:(k+1)*i]=aux
        b[j]=y[i]
        j=j+1
    end


    #Si(xi+1)=yi+1
    for i in 1:N-1
        aux=[x[i+1]^3, x[i+1]^2, x[i+1], 1]
        
        A[j,(k+1)*(i-1)+1:(k+1)*i]=aux
        b[j]=y[i+1]
        j=j+1
    end

    #S'i(xi+1)=S'i+1(xi+1)
    for i in 1:N-2
        aux=zeros(4*(N-1))
        aux[(k+1)*(i-1)+1:(k+1)*(i+1)]=[3*x[i+1]^2, 2*x[i+1], 1, 0, -3*x[i+1]^2, -2*x[i+1], -1, 0]      
        A[j,:]=aux
        b[j]=0
        j=j+1
    end

    #S''i(xi+1)=S''i+1(xi+1)
    for i in 1:N-2
        aux=[6*x[i+1], 2, 0, 0, -6*x[i+1], -2, 0,0]
        
        A[j,:]=aux
        b[j]=0
        j=j+1
    end

    #endpoints constraint
    aux=[6*x[1], 2, 0, 0, 0, 0, 0,0]
    A[j,:]=aux
    b[j]=0
    j=j+1
    aux=[0, 0, 0, 0, 6*x[N], 2, 0,0]
    A[j,:]=aux
    b[j]=0
     
    prob = LinearProblem(A, b)
    sol = solve(prob)

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

    println(cubic_f(x1,sol[(k+1)*(loc-1)+1:(k+1)*loc]))

    return sol[(k+1)*(loc-1)+1:(k+1)*loc]"""
end

function splinefour(x,y,x1)

    N=length(x)

    A=zeros(4*(N-1),4*(N-1))
    b=zeros(4*(N-1))
    k=4

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
    """for i in 1:N-2
        aux=[4*x[i+1]^3, 3*x[i+1]^2, 2*x[i+1], 1, 0,-4*x[i+1]^3, -3*x[i+1]^2, -2*x[i+1], -1,0]
        
        A[j,:]=aux
        b[j]=0
        j=j+1
    end

    #S''i(xi+1)=S''i+1(xi+1)
    for i in 1:N-2
        aux=[12*x[i+1]^2, 6*x[i+1], 2, 0, 0, -12*x[i+1]^2, -6*x[i+1], -2, 0, 0]
        
        A[j,:]=aux
        b[j]=0
        j=j+1
    end

    #endpoints constraint
    aux=[12*x[1]^3, 6*x[1], 2, 0, 0, 0, 0, 0, 0,0]
    A[j,:]=aux
    b[j]=0
    j=j+1
    aux=[0, 0, 0, 0, 0, 12*x[N]^3, 6*x[N], 2, 0,0]
    A[j,:]=aux
    b[j]=0
     
    prob = LinearProblem(A, b)
    sol = solve(prob)

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


x=[0 1 2 3]
y=[1 3 2 1]
using LinearSolve

x1=[0 1 2 3]
myspline(x,y,1.5)
#splinefour(x,y,1.5)



#evaluate(sol,x)