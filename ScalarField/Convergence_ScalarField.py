

import pandas as pd
import numpy as np
import os
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import style

import scipy
import numpy as np
from scipy import interpolate
from io import StringIO


# # Reading data for 3 resolutions

# Muninn format

# In[ ]:


datatime=[]
datagrid=[]
datam=[]
databeta=[]
datapsi=[]
dataderpsi=[]
vars=["m", "beta", "psi", "derpsi"]
res=[1,2,3]
for resolution in res:


    time=[]
    grid=[]
    m=[]
    beta=[]
    psi=[]
    derpsi=[]
    k=1

    for var in vars:
        #dir = "/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/StoredDATA/28-06/N200,dt=2.5e-5/res{}/{}.txt".format(resolution,var)
        #dir = "/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/DATA/muninnDATA/res{}/{}.txt".format(resolution,var)
        #dir = "/home/rita13santos/Desktop/29sepfieldplot/sub/{}.txt".format(var)
        #dir = "/home/rita13santos/Desktop/hojeconv/2ndeven/2ndeven/N200/res{}/{}.txt".format(resolution,var)
        dir = "/home/ritapsantos/data/ritapsantos/SavedDATA/conv2ndeven/res{}/{}.txt".format(resolution,var)
    
    
        print(dir)

        with open(dir) as f:
            for line in f:

                index = line.find("Time")
                if index==1:
                    if k==1:
                        time.append(float(line[index+7:len(line)-1]))
                    auxgrid=[]
                    auxdata=[]
                

                else:
                    a=line.split()
                    if a!=[]:
                        if k==1:
                            auxgrid.append(float(a[0]))
                        auxdata.append(float(a[1]))
                    elif a==[]:
                        grid.append(auxgrid)
                        if k==1:
                            m.append(auxdata)
                        elif k==2:
                            beta.append(auxdata)
                        elif k==3:
                            psi.append(auxdata)
                        elif k==4:
                            derpsi.append(auxdata)
                    
        k=k+1


    datatime.append(time)
    datagrid.append(grid)
    datam.append(m)
    databeta.append(beta)
    datapsi.append(psi)
    dataderpsi.append(derpsi)
#count=len(datatime[2])+1
#count


# In[ ]:


len(datagrid[0][0])


# In[ ]:


print(len(datam[0]))
print(len(databeta[0]))
print(len(datapsi[0]))
print(len(dataderpsi[0]))


# ####

# In[103]:


print(datam[0][0][1]) # datam indexes give res, time then gridpoint


# In[104]:


L=len(datam[0][0])-6 # grid length without the ghostpoints
dx=datagrid[0][0][1]-datagrid[0][0][0]


# In[105]:


#plt_x1 = np.linspace(0, 1, L)
#plt_x2 = np.linspace(0, 1, 2*L-1)
#plt_x3 = np.linspace(0, 1, 4*L-3)
t=0
plt_x1 = datagrid[0][t][3:len(datagrid[0][0])-3]
plt_x2 = datagrid[1][t*2][3:len(datagrid[1][0])-3]
plt_x3 = datagrid[2][t*3][3:len(datagrid[2][0])-3]

plt.plot(plt_x1,datam[0][0][3:len(datam[0][0])-3])
#plt.plot(plt_x2,datam[1][0][3:len(datam[1][0])-3])
#plt.plot(plt_x3,datam[2][0][3:len(datam[2][0])-3])


# # Styling plots

# In[106]:


plt.rcParams.update({
    'font.size': 12,
    'legend.fontsize':12,
    'xtick.labelsize': 'large',
    'xtick.color': 'black',
    'ytick.labelsize': 'large',
    'ytick.color': 'black'})


# # Plotting data with resolutions 1 and 2 and differences 

# In[107]:


# for given t
t1=0 #last timestep
t2=2*t1
t3=4*t1
auxm = []
auxbeta = []
auxpsi = []
auxderpsi = []
auxgrid = []

for i in range(len(datam[1][t2])):#iterate on the grid with higher resolution
    if ((i>2) and (i < (len(datam[1][t2])-3)) and ((i%2)!=0)): #ignoring ghost points
        auxm.append(datam[1][t2][i])
        auxbeta.append(databeta[1][t2][i])
        auxpsi.append(datapsi[1][t2][i])
        auxderpsi.append(dataderpsi[1][t2][i])
        auxgrid.append(datagrid[1][t2][i])
        #print(datam[1][0][i])

        
# array for all diff_m_2_1 diff_beta_2_1 diff_psi_2_1 diff_derpsi_2_1
diffs_2_1 = []

#subtracting f in both resolutions & ignoring ghost points
diffs_2_1.append(np.subtract(datam[0][t1][3:len(datam[0][t1])-3],auxm))
diffs_2_1.append(np.subtract(databeta[0][t1][3:len(databeta[0][t1])-3],auxbeta))
diffs_2_1.append(np.subtract(datapsi[0][t1][3:len(datapsi[0][t1])-3],auxpsi))
diffs_2_1.append(np.subtract(dataderpsi[0][t1][3:len(dataderpsi[0][t1])-3],auxderpsi))
diffs_2_1.append(np.subtract(datagrid[0][t1][3:len(datagrid[0][t1])-3],auxgrid))

with plt.style.context('ggplot'):
    plt.plot(plt_x1, diffs_2_1[0], label = 'res2-res1 of m')

plt.legend()
plt.xlabel('x')
plt.ylabel('m2(x)-m1(x)')
#plt.show()


# # Pointwise convergence tests: plots for thesis project report 

# In[108]:


# for given t

aux2m = []
aux2beta = []
aux2psi = []
aux2derpsi = []
aux2grid = []
for i in range(len(datam[2][t3])):#iterate on the grid with higher resolution
    if ((i>2) and (i < (len(datam[2][t3])-3)) and (((i+1)%4)==0)): #ignoring ghost points
        aux2m.append(datam[2][t3][i])
        aux2beta.append(databeta[2][t3][i])
        aux2psi.append(datapsi[2][t3][i])
        aux2derpsi.append(dataderpsi[2][t3][i])
        aux2grid.append(datagrid[2][t3][i])
        


#subtracting m beta psi and derpsi in both resolutions & ignoring ghost points
diffs_4_2 =[]
diffs_4_2.append(np.subtract(auxm,aux2m))
diffs_4_2.append(np.subtract(auxbeta,aux2beta))
diffs_4_2.append(np.subtract(auxpsi,aux2psi))
diffs_4_2.append(np.subtract(auxderpsi,aux2derpsi))
diffs_4_2.append(np.subtract(auxgrid,aux2grid))

with plt.style.context('ggplot'):
    
    fig, ax = plt.subplots(2, 2,figsize=(15, 10))

    ax[0][0].plot(plt_x1, (diffs_2_1[0]), label='Res1-res2 of m')
    ax[0][0].plot(plt_x1, (diffs_4_2[0])*4, label='Res4-res2 * 16 of m')
    ax[0][0].legend()
    
    ax[0][1].plot(plt_x1, (diffs_2_1[1]), label='Res1-res2 of beta')
    ax[0][1].plot(plt_x1, (diffs_4_2[1])*4, label='Res4-res2 *16 of beta')
    #ax[0][1].set_xlim([0, 0.05])
    #ax[0][1].set_ylim([0, 1*10**(-11)])
    ax[0][1].legend()
    
    ax[1][0].plot(plt_x1, (diffs_2_1[2]), label='Res1-res2 of psi')
    ax[1][0].plot(plt_x1, (diffs_4_2[2])*4, label='Res4-res2 * 16 of psi')
    ax[1][0].set_xlim([0, 1])
    ax[1][0].legend()
    
    ax[1][1].plot(plt_x1, (diffs_2_1[3]), label='Res1-res2 of derpsi')
    ax[1][1].plot(plt_x1, (diffs_4_2[3])*4, label='Res4-res2 * 16 of derpsi')
    #ax[1][1].set_xlim([0, 0.2])
    #ax[1][1].set_ylim([0, 0.5*10**(-9)])
    ax[1][1].legend()
    
    """ax[1][1].plot(plt_x1, abs(diffs_2_1[4]), label='Res1-res2 of spline derivative')
    ax[1][1].plot(plt_x1, abs(diffs_4_2[4])*16, label='Res4-res2 * 16 of spline derivative')
    #ax[1][1].set_xlim([0, 0.2])
    #ax[1][1].set_ylim([0, 0.5*10**(-9)])
    ax[1][1].legend()"""

plt.legend()
#plt.show()


# # L2 norm convergence tests

# In[109]:


#for the first time step

norm_low_med = []
norm_med_high = []

for j in range(5):
    a = 0
    b = 0
    for i in diffs_2_1[j]:
        a += i**2
    for i in diffs_4_2[j]:
        b += i**2
    norm_low_med.append(np.sqrt(a))
    norm_med_high.append(np.sqrt(b))

print("m convergence factor:")
print(math.log2(norm_low_med[0]/norm_med_high[0]))
print("")

print("beta convergence factor:")
print(math.log2(norm_low_med[1]/norm_med_high[1]))
print("")

print("psi convergence factor:")
print(math.log2(norm_low_med[2]/norm_med_high[2]))
print("")

print("derpsi convergence factor:")
print(math.log2(norm_low_med[3]/norm_med_high[3]))
print("")


# # Q(t)

# In[110]:


final_t=len(datatime[2])


# In[111]:


"""m200=Q_m
beta200=Q_beta
psi200=Q_psi
derpsi200=Q_derpsi
x200=plt_x = np.linspace(0, 5.0, len(m200))
Q200=Q""";


# In[112]:


rsquarednorm=True


# In[118]:


"""from cycler import cycler
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))"""

Q_m = []
Q_beta = []
Q_psi = []
Q_derpsi = []
Q= []


for t in range(0,int(final_t/4)):
    # for given t
    t1=t #last timestep
    t2=2*t1
    t3=4*t1
    auxm = []
    auxbeta = []
    auxpsi = []
    auxderpsi = []

    for i in range(len(datam[1][t2])):#iterate on the grid with higher resolution
        if ((i>2) and (i < (len(datam[1][t2])-3)) and ((i%2)!=0)): #ignoring ghost points
            auxm.append(datam[1][t2][i])
            auxbeta.append(databeta[1][t2][i])
            auxpsi.append(datapsi[1][t2][i])
            auxderpsi.append(dataderpsi[1][t2][i])
            #print(datam[1][0][i])


    # array for all diff_m_2_1 diff_beta_2_1 diff_psi_2_1 diff_derpsi_2_1
    diffs_2_1 = []

    #subtracting f in both resolutions & ignoring ghost points
    diffs_2_1.append(np.subtract(datam[0][t1][3:len(datam[0][t1])-3],auxm))
    diffs_2_1.append(np.subtract(databeta[0][t1][3:len(databeta[0][t1])-3],auxbeta))
    diffs_2_1.append(np.subtract(datapsi[0][t1][3:len(datapsi[0][t1])-3],auxpsi))
    diffs_2_1.append(np.subtract(dataderpsi[0][t1][3:len(dataderpsi[0][t1])-3],auxderpsi))

    # for given t

    aux2m = []
    aux2beta = []
    aux2psi = []
    aux2derpsi = []
    for i in range(len(datam[2][t3])):#iterate on the grid with higher resolution
        if ((i>2) and (i < (len(datam[2][t3])-3)) and (((i+1)%4)==0)): #ignoring ghost points
            aux2m.append(datam[2][t3][i])
            aux2beta.append(databeta[2][t3][i])
            aux2psi.append(datapsi[2][t3][i])
            aux2derpsi.append(dataderpsi[2][t3][i])



    #subtracting m beta psi and derpsi in both resolutions & ignoring ghost points
    diffs_4_2 =[]
    diffs_4_2.append(np.subtract(auxm,aux2m))
    diffs_4_2.append(np.subtract(auxbeta,aux2beta))
    diffs_4_2.append(np.subtract(auxpsi,aux2psi))
    diffs_4_2.append(np.subtract(auxderpsi,aux2derpsi))


    ####

    norm_low_med = []
    norm_med_high = []
    single_norm_low_med = []
    single_norm_med_high = []
    s_a=0
    s_b=0

    for j in range(4):
        a = 0
        b = 0
        for i in diffs_2_1[j]:
            a += i**2
            s_a += i**2
        for i in diffs_4_2[j]:
            b += i**2
            s_b += i**2
        norm_low_med.append(np.sqrt(a))
        norm_med_high.append(np.sqrt(b))

    single_norm_low_med=np.sqrt(s_a)
    single_norm_med_high=np.sqrt(s_b)
    
    Q_m.append(math.log2(norm_low_med[0]/norm_med_high[0]))
    Q_beta.append(math.log2(norm_low_med[1]/norm_med_high[1]))
    Q_psi.append(math.log2(norm_low_med[2]/norm_med_high[2]))
    Q_derpsi.append(math.log2(norm_low_med[3]/norm_med_high[3]))
    Q.append(math.log2(single_norm_low_med/single_norm_med_high))
    ####
    
plt_x = np.linspace(0, 5.0, len(Q_m))
#plt_x = np.linspace(0, int(final_t/4)*dx*step, int(final_t/4))
with plt.style.context('ggplot'):
    
    fig, ax = plt.subplots(1, 2,figsize=(10, 5))
    
    ax[0].plot(plt_x, Q_m, label = '$Q_m$',color='#1f77b4')
    ax[0].plot(plt_x, Q_beta, label = '$Q_{\\beta}$',color='#ff7f0e')
    ax[0].plot(plt_x, Q_psi, label = '$Q_{\psi}$',color='#2ca02c')
    ax[0].plot(plt_x, Q_derpsi, label = '$Q_{\psi,r}$',color='#d62728')
    
    
    ax[0].plot(x200, m200, linestyle='dotted',alpha=0.7,color='#1f77b4')
    ax[0].plot(x200, beta200, linestyle='dotted',alpha=0.7,color='#ff7f0e')
    ax[0].plot(x200, psi200, linestyle='dotted',alpha=0.7,color='#2ca02c')
    ax[0].plot(x200, derpsi200, linestyle='dotted',alpha=0.7,color='#d62728')
    
    ax[0].set_ylim([0,8])
    ax[0].set_xlim([0,2.5])
    ax[0].legend(loc ="upper left")
    ax[0].set_xlabel('Proper time u')
    ax[0].set_ylabel('Q(u)')
    ax[0].set_title('Individual norms')
    
    
    ax[1].plot(plt_x, Q, label = '$Q$',color='#1f77b4')
    ax[1].plot(x200, Q200, linestyle='dotted',alpha=0.7,color='#1f77b4')
    #ax[1].plot(x150, Q150, linestyle='dotted',alpha=0.7,color='#1f77b4')
    ax[1].set_ylim([0,8])
    ax[1].set_xlim([0,2.5])
    ax[1].legend(loc ="upper left")
    ax[1].set_title('Single norm')
    #plt.xticks(rotation=45)
    ax[1].set_xlabel('Proper time u')
    ax[1].set_ylabel('Q(t)')

plt.savefig("Qt.pdf", format="pdf", bbox_inches="tight")
#plt.show()


# In[114]:

