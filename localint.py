# Python code accompanying 'How local interactions impact the dynamics of an epidemic' by Wren & Best
# This code will reproduce figure 3a from that manuscript
# Code developed by Alex Best and Lydia Wren in 2020

# Import the necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clt
from scipy.integrate import odeint
from random import *

# Constants
N=25           # Grid aide length
tot=N*N        # Total Grid size
GAMMA=1/14     # Recovery rate
R0=2           # R0
BETA=R0*GAMMA/tot   # Transmission coefficient
REPS=100       # Simulation replicates 
NCURVES=20     # Number of curves to use for sampling
NSAMPLES=50    # Number of samples
I0=5           # Initial number of infecteds

bigstore=[]
bigstore2=[]
bigstored=[]

# Functions

# Function to count noumbers in each compartment
def count_type(type):
    return np.count_nonzero(grid[:,:] == type)

# Function to find nearest-neighbours
def local_inf(site_x,site_y):
    local_inf=0
    if site_x==0:
        downx=N-1
    else:
        downx=site_x-1
    if site_y==0:
        downy=N-1
    else:
        downy=site_y-1
    if site_x==N-1:
        upx=0
    else:
        upx=site_x+1
    if site_y==N-1:
        upy=0
    else:
        upy=site_y+1
        
    if grid[downx,site_y]==1:
        local_inf=local_inf+1
    if grid[upx,site_y]==1:
        local_inf=local_inf+1
    if grid[site_x,downy]==1:
        local_inf=local_inf+1
    if grid[site_x,upy]==1:
        local_inf=local_inf+1
    return local_inf

# Function to check the current scale
def findscale():
    S=count_type(0)
    I=count_type(1)
    localinf=0
    for i in range(0,N):
        for j in range(0,N):
            if grid[i,j]==0:
                localinf=local_inf(i,j)+localinf
    if S>0:
        QIS=localinf/S*tot/4
    else:
        QIS=0
    #Set relative parameter values
    scale=GAMMA*I+BETA*S*(L*QIS+(1-L)*I)  
    return scale

# MAIN ROUTINE

for ll in range (0,11):     # First run sets L=LLow, Second run has control with Llow and Lhgih
    L=ll*0.1
    print(L)
    peaks=[] # For storing peak no. of infecteds
    tots=[] # For storing total no. infected
    days=[]
    for reps in range(0,REPS):
        # Set initial conditions
        init_inf=I0
        grid=np.zeros((N,N))
        for i in range(0,init_inf):
            grid[randint(0,N-1),randint(0,N-1)]=1

        tsteps=[0]
        infecteds=[count_type(1)]
        current_t=0

        # Main run
        while current_t<300:
            # Find tau-leap         
            scale=findscale()
            dt = -np.log(random()) / scale
            current_t=tsteps[-1]

            # Create randomised list of sites to check
            findx=[i for i in range(N)]
            findy=[i for i in range(N)]
            shuffle(findx)
            shuffle(findy) 
            flagged=0   # Used to break out of 2nd loop

            #Find event
            if random()<GAMMA*infecteds[-1]/scale: #Event is recovery   
                for tryx in findx:
                    if flagged==1:
                        break
                    for tryy in findy:
                        if grid[tryx,tryy]==1:
                            grid[tryx,tryy]=2
                            flagged=1
                            break
            else: #Event is transmission
                if random()>L: #Transmission is global
                    for tryx in findx:
                        if flagged==1:
                            break
                        for tryy in findy:
                            if grid[tryx,tryy]==0:                       
                                grid[tryx,tryy]=1
                                flagged=1
                                break
                else: # Transmission is local
                    for tryx in findx:
                        if flagged==1:
                            break
                        for tryy in findy:
                            if grid[tryx,tryy]==0: 
                                 if local_inf(tryx,tryy)>0:
                                    grid[tryx,tryy]=1
                                    flagged=1
                                    break

            # Update time and infection lists
            tsteps.append(dt+current_t)
            infecteds.append(count_type(1))
            if infecteds[-1]==0:
                break
        
        peaks.append(max(infecteds)/tot)
        day=np.where(np.array(infecteds)==max(infecteds))
        days.append(tsteps[day[0][0]])
        tots.append((count_type(1)+count_type(2))/tot)
    
    print("The median peak is", np.median(peaks))
    print("The median total is", np.median(tots))
    
    bigstore.append(peaks)
    bigstore2.append(tots)
    bigstored.append(days)

yy=[j*0.1 for j in range (0,11)]
xx = [ '%.1f' % elem for elem in yy ]
plt.rcParams.update({'font.size': 16})

fig, ax = plt.subplots()
data=[bigstore[:][0],bigstore[:][1],bigstore[:][2],bigstore[:][3],bigstore[:][4],bigstore[:][5],bigstore[:][6],bigstore[:][7],bigstore[:][8],bigstore[:][9],bigstore[:][10]]
plt.boxplot(data)
ax.set(xticklabels=xx)
plt.xlabel('$L$')
plt.ylabel('Peak infected')
plt.ylim(0,0.3)
plt.tight_layout()
plt.savefig('peaksbox.png')

fig, ax = plt.subplots()
data=[bigstore2[:][0],bigstore2[:][1],bigstore2[:][2],bigstore2[:][3],bigstore2[:][4],bigstore2[:][5],bigstore2[:][6],bigstore2[:][7],bigstore2[:][8],bigstore2[:][9],bigstore2[:][10]]
plt.boxplot(data)
ax.set(xticklabels=xx)
plt.xlabel('$L$')
plt.ylabel('Total infected')
plt.ylim(0,1)
plt.tight_layout()
plt.savefig('totsbox.png')

fig, ax = plt.subplots()
data=[bigstored[:][0],bigstored[:][1],bigstored[:][2],bigstored[:][3],bigstored[:][4],bigstored[:][5],bigstored[:][6],bigstored[:][7],bigstored[:][8],bigstored[:][9],bigstored[:][10]]
plt.boxplot(data)
ax.set(xticklabels=xx)
plt.xlabel('$L$')
plt.ylabel('Day of Peak')
plt.ylim(0,300)
plt.tight_layout()
plt.savefig('daysbox.png')

np.savetxt("daystore2.txt",bigstored,fmt="%.3f")
np.savetxt("totstore2.txt",bigstore2,fmt="%.3f")
np.savetxt("peakstore2.txt",bigstore,fmt="%.3f")