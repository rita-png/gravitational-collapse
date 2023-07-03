#!/usr/bin/env python
# coding: utf-8

# In[43]:


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
        dir = "/home/rita13santos/Desktop/MSc Thesis/Git/ScalarField/DATA/muninnDATA/res{}/{}.txt".format(resolution,var)
        #dir = "/home/rita13santos/Desktop/muninnDATA/res{}/{}.txt".format(resolution,var)
    
    
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


print(len(datam[0]))
print(len(databeta[0]))
print(len(datapsi[0]))
print(len(dataderpsi[0]))




print(datam[0][0][1]) # datam indexes give res, time then gridpoint


# In[35]:


L=len(datam[0][0])-6 # grid length without the ghostpoints
dx=datagrid[0][0][1]-datagrid[0][0][0]


# In[36]:


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

# In[37]:


plt.rcParams.update({'font.size': 12})


# # Plotting data with resolutions 1 and 2 and differences 

# In[38]:


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

with plt.style.context('bmh'):
    plt.plot(plt_x1, diffs_2_1[0], label = 'res2-res1 of m')

plt.legend()
plt.xlabel('x')
plt.ylabel('m2(x)-m1(x)')
plt.show()


# # Pointwise convergence tests: plots for thesis project report 

# In[39]:


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

"""with plt.style.context('bmh'):
    
    fig, ax = plt.subplots(2, 2,figsize=(15, 10))

    ax[0][0].plot(plt_x1, (diffs_2_1[0]), label='Res1-res2 of m')
    ax[0][0].plot(plt_x1, (diffs_4_2[0])*16, label='Res4-res2 * 16 of m')
    ax[0][0].legend()
    
    ax[0][1].plot(plt_x1, (diffs_2_1[1]), label='Res1-res2 of beta')
    ax[0][1].plot(plt_x1, (diffs_4_2[1])*16, label='Res4-res2 *16 of beta')
    #ax[0][1].set_xlim([0, 0.05])
    #ax[0][1].set_ylim([0, 1*10**(-11)])
    ax[0][1].legend()
    
    ax[1][0].plot(plt_x1, (diffs_2_1[2]), label='Res1-res2 of psi')
    ax[1][0].plot(plt_x1, (diffs_4_2[2])*16, label='Res4-res2 * 16 of psi')
    ax[1][0].set_xlim([0, 1])
    ax[1][0].legend()
    
    ax[1][1].plot(plt_x1, (diffs_2_1[3]), label='Res1-res2 of derpsi')
    ax[1][1].plot(plt_x1, (diffs_4_2[3])*16, label='Res4-res2 * 16 of derpsi')
    #ax[1][1].set_xlim([0, 0.2])
    #ax[1][1].set_ylim([0, 0.5*10**(-9)])
    ax[1][1].legend()
    

plt.legend()
plt.show()"""


# # L2 norm convergence tests

# In[40]:


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

# In[41]:


final_t=len(datatime[2])


# In[42]:


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
    
plt_x = np.linspace(0, 3.0, len(Q_m))
#plt_x = np.linspace(0, int(final_t/4)*dx*step, int(final_t/4))
with plt.style.context('bmh'):
    
    fig, ax = plt.subplots(1, 2,figsize=(10, 5))
    
    ax[0].plot(plt_x, Q_m, label = '$Q_m$')
    ax[0].plot(plt_x, Q_beta, label = '$Q_{\\beta}$',alpha=0.5)
    ax[0].plot(plt_x, Q_psi, label = '$Q_{\psi}$')
    ax[0].plot(plt_x, Q_derpsi, label = '$Q_{\psi,x}$',alpha=0.5)
    ax[0].set_ylim([0,8])
    ax[0].set_xlim([0,3.0])
    ax[0].legend(loc ="upper left")
    ax[0].set_xlabel('time')
    ax[0].set_ylabel('Q(t)')
    ax[0].set_title('Individual norms')
    
    
    ax[1].plot(plt_x, Q, label = '$Q$',alpha=0.5)
    ax[1].set_ylim([0,8])
    ax[1].set_xlim([0,3.0])
    ax[1].legend(loc ="upper left")
    ax[1].set_title('Single norm')
    plt.xticks(rotation=45)
    ax[1].set_xlabel('time')
    ax[1].set_ylabel('Q(t)')

plt.show()


# In[13]:


Q_m


# # Pointwise convergence through time

# In[14]:


time_frame=range(0,int(final_t/4))
pw_m_21 = []
pw_m_42 = []
pw_beta_21 = []
pw_beta_42 = []
pw_psi_21 = []
pw_psi_42 = []
pw_derpsi_21 = []
pw_derpsi_42 = []

for i in time_frame:
    # for given t
    t1=i #last timestep
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

    pw_m_21.append(abs(diffs_2_1[0]))
    pw_m_42.append(abs(diffs_4_2[0])*16)
    pw_beta_21.append(abs(diffs_2_1[1]))
    pw_beta_42.append(abs(diffs_4_2[1])*16)
    pw_psi_21.append(abs(diffs_2_1[2]))
    pw_psi_42.append(abs(diffs_4_2[2])*16)
    pw_derpsi_21.append(abs(diffs_2_1[3]))
    pw_derpsi_42.append(abs(diffs_4_2[3])*16)


# In[16]:


step=5


time_frame=range(0,int(final_t/4))
pw_m_21 = []
pw_m_42 = []
pw_beta_21 = []
pw_beta_42 = []
pw_psi_21 = []
pw_psi_42 = []
pw_derpsi_21 = []
pw_derpsi_42 = []

for i in time_frame:
    # for given t
    t1=i #last timestep
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

    pw_m_21.append(abs(diffs_2_1[0]))
    pw_m_42.append(abs(diffs_4_2[0])*16)
    pw_beta_21.append(abs(diffs_2_1[1]))
    pw_beta_42.append(abs(diffs_4_2[1])*16)
    pw_psi_21.append(abs(diffs_2_1[2]))
    pw_psi_42.append(abs(diffs_4_2[2])*16)
    pw_derpsi_21.append(abs(diffs_2_1[3]))
    pw_derpsi_42.append(abs(diffs_4_2[3])*16)


# In[16]:


step=5


# In[18]:


from matplotlib import pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots(2,2,figsize=(15, 10))
for j in range(0,2):
    for i in range(0,2):
        ax[i][j].set_xlim(0, 1)
        ax[i][j].grid()
#2th
ax[0][0].set_ylim(0, 2*10**(-9))
ax[0][1].set_ylim(0, 10**(-8))
ax[1][0].set_ylim(0, 10**(-7))
ax[1][1].set_ylim(0, 5*10**(-7))
#4th
#ax[0][0].set_ylim(0, 2*10**(-11))
#ax[0][1].set_ylim(0, 10**(-8))
#ax[1][0].set_ylim(0, 10**(-9))
#ax[1][1].set_ylim(0, 10**(-9))

line1, = ax[0][0].plot([], [], lw = 3)
line2, = ax[0][0].plot([], [], lw = 3)

time_text = ax[0][0].text(0.5, 0.9, '', transform=ax[0][0].transAxes)

line3, = ax[0][1].plot([], [], lw = 3)
line4, = ax[0][1].plot([], [], lw = 3)

line5, = ax[1][0].plot([], [], lw = 3)
line6, = ax[1][0].plot([], [], lw = 3)

line7, = ax[1][1].plot([], [], lw = 3)
line8, = ax[1][1].plot([], [], lw = 3)


def init():
    line1.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    return line1,line2,time_text

def animate(i):
    # m
    x = plt_x1
    y = pw_m_21[i]
    y2 = pw_m_42[i]
    line1.set_data(x, y)
    line2.set_data(x, y2)
    line1.set_label('m21')
    line2.set_label('m42 *4')
    ax[0][0].legend()
    
    # beta
    x = plt_x1
    y = pw_beta_21[i]
    y2 = pw_beta_42[i]
    line3.set_data(x, y)
    line4.set_data(x, y2)
    line3.set_label('beta21')
    line4.set_label('beta42 *4')
    ax[0][1].legend()
    
    # psi
    x = plt_x1
    y = pw_psi_21[i]
    y2 = pw_psi_42[i]
    line5.set_data(x, y)
    line6.set_data(x, y2)
    line5.set_label('psi21')
    line6.set_label('psi42 *4')
    ax[1][0].legend()
    
    # der psi
    x = plt_x1
    y = pw_derpsi_21[i]
    y2 = pw_derpsi_42[i]
    line7.set_data(x, y)
    line8.set_data(x, y2)
    line7.set_label('derpsi21')
    line8.set_label('derpsi42 *16')
    ax[1][1].legend()
    
    t=i*dx*step
    time_text.set_text('time = %.3f' % t)
    
    return line1,line2,line3,line4,line5,line6,line7,line8,time_text
 

anim = FuncAnimation(fig, animate,
                    init_func = init,
                    frames = len(time_frame),
                    interval = 200,
                    blit = False,repeat=False)

anim.save('pointwise.gif',
          writer = 'ffmpeg', fps = 2*4)
