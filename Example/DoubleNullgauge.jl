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


#Building initial data with a  Runge-Kutta integrator for the constraint

function rungekutta4(f,y0,T)
    n = length(T)
    y = zeros((n, length(y0)))
    y[1,:] = y0;
    for i in 1:n-1
        h = T[2] - T[1]
        k1 = f(y[i,:], T[i])
        k2 = f(y[i,:] + k1 * h/2, T[i] + h/2)
        k3 = f(y[i,:] + k2 * h/2, T[i] + h/2)
        k4 = f(y[i,:] + k3 * h, T[i] + h)
        y[i+1,:] = y[i,:] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    end
    return y
end


# Runge Kutta integrator used for the method of lines

function rungekutta4molstep(f,y0,T,w::Int64,ex)
    y = zeros((length(R),length(y0[1,:])));
    y[:,:] = y0[:,:]
        h = dt; # only for equally spaced grid in time and space, otherwise (T[w+1] - T[w])
        k1 = f(y[:,:], T[w])
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

#ghosts

function ghost(y)
    L=length(y[:,1])
    y[1,1]=y[5,1]; #delta is even
    y[2,1]=y[4,1];
    y[1,2]=y[5,3]; #sdelta
    y[2,2]=y[4,3];
    y[1,3]=y[5,2]; #sbardelta
    y[2,3]=y[4,2];
    y[1,4]=-y[5,4]; #Rc is odd
    y[2,4]=-y[4,4];
    y[1,5]=-y[5,6]; #sRc 
    y[2,5]=-y[4,6];
    y[1,6]=-y[5,5]; #sbarRc
    y[2,6]=-y[4,5];
    y[1,7]=y[5,7]; #Pi is even
    y[2,7]=y[4,7];
    y[1,8]=-y[5,8]; # Phi is odd
    y[2,8]=-y[4,8];
    y[1,9]=y[5,9]; #phi is even
    y[2,9]=y[4,9];    
    #y[2,:]=2*y[3,:]-y[4,:];
    #y[1,:]=2*y[2,:]-y[3,:];
    y[L-1,:]=2*y[L-2,:]-y[L-3,:]; #extrapolation in the outer boundary
    y[L,:]=2*y[L-1,:]-y[L-2,:];
    return y
end

#fourth order  dissipation

function dissipation4(y,i)
        delta4=(y[i+2,:]-4*y[i+1,:]+6*y[i,:]-4*y[i-1,:]+y[i-2,:]);
    return (-1)^4*epsilon*1/(dx)*delta4
end

DNconstraint(y,R1)=[(1.0-y[1]^2.0)/(2.0*y[2])-4.0*pi*y[2]*scalar_spaceder(R1)^2.0 y[1]]; 

function DNconstraintRHS(y,R1)
        if R1<10^(-7)
                return [0.0 1.0]'
        else
                return DNconstraint(y,R1)'
        end
end

# Discretization of derivatives
Der(y,i,k)=(y[i+1,k]-y[i-1,k])/(R[i+1]-R[i-1]);
DDer(y,i,k)=(y[i+1,k]-2.0*y[i,k]+y[i-1,k])/(R[i+1]-R[i-1])^2.0;

# RHSs for the bulk equations

sRcRHS(y,i)=exp(-y[i,1])*(-4.0*pi*y[i,4]*(y[i,7]+y[i,8])^2.0)-Der(y,i,5);       #sRc
sbarRcRHS(y,i)=exp(-y[i,1])*(-4.0*pi*y[i,4]*(y[i,7]-y[i,8])^2.0)+Der(y,i,6);     	#sbarRc
sdRHS(y,i)=(2.0+2.0*exp(y[i,1])*y[i,6]*y[i,5])/(y[i,4]^2.0)-exp(y[i,1])*y[i,3]*y[i,2]-8.0*exp(-y[i,1])*pi*y[i,4]^2.0*(y[i,7]^2.0-y[i,8]^2.0)+Der(y,i,2);		#sd
sbardRHS(y,i)=(2.0+2.0*exp(y[i,1])*y[i,6]*y[i,5])/(y[i,4]^2.0)-exp(y[i,1])*y[i,3]*y[i,2]-8.0*exp(-y[i,1])*pi*y[i,4]^2.0*(y[i,7]^2.0-y[i,8]^2.0)-Der(y,i,3);		#sbard
ScalarRHSDN(y,i)=(-exp(y[i,1])*y[i,7]*(y[i,6]+y[i,5]))/y[i,4]+3.0*Der(y,i,4)/(y[i+1,4]^3.0-y[i-1,4]^3.0)*(y[i+1,8]*y[i+1,4]^2.0-y[i-1,8]*y[i-1,4]^2.0);  	#Pi 

OriScalarRHSDN(y,i)=-exp(y[i,1])*y[i,7]*EvansScalar(y,i)+3.0*Der(y,i,4)/(y[i+1,4]^3.0-y[i-1,4]^3.0)*(y[i+1,8]*y[i+1,4]^2.0-y[i-1,8]*y[i-1,4]^2.0);

TaylorMms(y,i)=exp(y[i,1])*(2.0*Der(y,i,5)*Der(y,i,6)+y[i,5]*(2.0*Der(y,i,1)*Der(y,i,6)+DDer(y,i,6))+y[i,6]*(2.0*Der(y,i,1)*Der(y,i,5)+y[i,5]*(Der(y,i,1)^2.0+DDer(y,i,1))+DDer(y,i,5)))/(Der(y,i,4)^2.0+y[i,4]*DDer(y,i,4));

EvansScalar(y,i)=3.0*((y[i+1,5]+y[i+1,6])*R[i+1]^2.0-(y[i-1,5]+y[i-1,6])*R[i-1]^2.0)/(R[i+1]^3.0-R[i-1]^3.0)-Der(y,i,5)-Der(y,i,6);

EvansMms(y,i)=(2.0+2.0*exp(y[i,1])*y[i,6]*y[i,5])*3.0*Der(y,i,4)*(R[i+1]-R[i-1])/(y[i+1,4]^3.0-y[i-1,4]^3.0);

function originDN(y,i)
    dy=zeros(length(y[1,:]));
    dy[1]=exp(y[i,1])*y[i,3]+Der(y,i,1)-dissipation4(y,i)[1]; #delta
    dy[2]=EvansMms(y,i)-exp(y[i,1])*y[i,3]*y[i,2]-8.0*exp(-y[i,1])*pi*y[i,4]^2.0*(y[i,7]^2.0-y[i,8]^2.0)+Der(y,i,2)-dissipation4(y,i)[2]; #sdelta
    dy[3]=EvansMms(y,i)-exp(y[i,1])*y[i,3]*y[i,2]-8.0*exp(-y[i,1])*pi*y[i,4]^2.0*(y[i,7]^2.0-y[i,8]^2.0)-Der(y,i,3)-dissipation4(y,i)[3]; 	#sbardelta
    dy[4]=exp(y[i,1])*y[i,6]+Der(y,i,4)-dissipation4(y,i)[4];	#Rc
    dy[5]=sRcRHS(y,i)-dissipation4(y,i)[5];	#sRc
    dy[6]=sbarRcRHS(y,i)-dissipation4(y,i)[6];	#sbarRc
    dy[7]=OriScalarRHSDN(y,i)-dissipation4(y,i)[7];  #Pi
    dy[8]=Der(y,i,7)-dissipation4(y,i)[8];	#Phi
    dy[9]=y[i,7]-dissipation4(y,i)[9];		#phi
    return dy
end


function bulkDN(y,i)
    dy=zeros(length(y[1,:]));
    dy[1]=exp(y[i,1])*y[i,3]+Der(y,i,1)-dissipation4(y,i)[1]; #delta
    dy[2]=sdRHS(y,i)-dissipation4(y,i)[2]; #sdelta
    dy[3]=sbardRHS(y,i)-dissipation4(y,i)[3];   #sbardelta
    dy[4]=exp(y[i,1])*y[i,6]+Der(y,i,4)-dissipation4(y,i)[4];   #Rc
    dy[5]=sRcRHS(y,i)-dissipation4(y,i)[5];     #sRc
    dy[6]=sbarRcRHS(y,i)-dissipation4(y,i)[6];  #sbarRc
    dy[7]=ScalarRHSDN(y,i)-dissipation4(y,i)[7];  #Pi
    dy[8]=Der(y,i,7)-dissipation4(y,i)[8];      #Phi
    dy[9]=y[i,7]-dissipation4(y,i)[9];          #phi
    return dy
end

function boundaryDN(y,i)
    dy=zeros(length(y[1,:]));
    dy[1]=exp(y[i,1])*y[i,3]+Der(y,i,1)-dissipation4(y,i)[1]; #delta
    dy[2]=sdRHS(y,i)-dissipation4(y,i)[2]; #sdelta
    dy[3]=sbardRHS(y,i)-dissipation4(y,i)[3];   #sbardelta
    dy[4]=exp(y[i,1])*y[i,6]+Der(y,i,4)-dissipation4(y,i)[4];   #Rc
    dy[5]=sRcRHS(y,i)-dissipation4(y,i)[5];     #sRc
    dy[6]=sbarRcRHS(y,i)-dissipation4(y,i)[6];  #sbarRc
    dy[7]=ScalarRHSDN(y,i)-dissipation4(y,i)[7];  #Pi
    dy[8]=Der(y,i,7)-dissipation4(y,i)[8];      #Phi
    dy[9]=y[i,7]-dissipation4(y,i)[9];          #phi
    return dy
end

# Defining the function in the RHS of the evolution equation system

function DNRHS(y,T)
    L=length(R)
    dy=zeros((L,length(y[1,:])));
        for i in 3:(L-2)
                dy[i,:]=originDN(y,i);
        end
    #dy[3,:]=originDN(y,3);
    #dy[L-2,:]=boundaryDN(y,L-2);
    return dy
end


function fixedbulkDN(y,i)
    #y=ghost(y);
    dy=zeros(length(y[1,:]));
    dy[1]=0; #exp(y[i,1])*y[i,3]+Der(y,i,1)-dissipation4(y,i)[1]; #delta
    dy[2]=0;  #sdRHS(y,i)-dissipation4(y,i)[2]; #sdelta
    dy[3]=0; #sbardRHS(y,i)-dissipation4(y,i)[3];   #sbardelta
    dy[4]=0; #exp(y[i,1])*y[i,6]+Der(y,i,4)-dissipation4(y,i)[4];   #Rc
    dy[5]=0; #sRcRHS(y,i)-dissipation4(y,i)[5];     #sRc
    dy[6]=0; #sbarRcRHS(y,i)-dissipation4(y,i)[6];  #sbarRc
    dy[7]=OriScalarRHSDN(y,i)-dissipation4(y,i)[7];  #Pi
    dy[8]=Der(y,i,7)-dissipation4(y,i)[8];      #Phi
    dy[9]=y[i,7]-dissipation4(y,i)[9];          #phi
    return dy
end

function fixedDNRHS(y,T)
    #dx=R[i+1]-R[i]
    L=length(R)
    dy=zeros((L,length(y[1,:])));
        for i in 3:(L-2)
                dy[i,:]=fixedbulkDN(y,i);
        end
    #dy[3,:]=originDN(y,L-2);
    #dy[L-2,:]=boundaryDN(y,L-2);
    return dy
end

