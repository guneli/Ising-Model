from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
import os.path
## change these parameters for a smaller (faster) simulation 
nt      = 88     #  number of temperature points
eqSteps = 3000  #  number of MC sweeps for equilibration
mcSteps = 20000     #  number of MC sweeps for calculation
L=[8,16,35]#System size
T = np.linspace(1.53, 3.28, nt); 
mc= np.linspace(100,mcSteps, 200); #only 100th configuration will be considered
configurations=[]
El=[];
Ml=[];
Cl=[];
Xl=[];
BCl=[];
def writeToFile(T,E,M,C,X,BC,L):
    fileName ="Properties_L_{}".format(L)
    completeName = os.path.join("../Ising/", fileName+".txt")

    if os.path.isfile(completeName):
        file= open(completeName,"a")
        file.write("{} \n {} \n {} \n {}  \n {}  \n {} \n".format(T,E,M,C,X,BC))
    else:
        file= open(completeName,"a")
        file.write("{} \n {} \n {} \n {}  \n {} \n {}  \n".format(T,E,M,C,X,BC))
    file.close()
def init_lattice(n):

    '''Create a nxn lattice with random spin'''
    
    lattice = np.random.choice([1, -1], size=(n, n))
    return lattice

def mcmove(lattice, beta,N):
    '''Monte Carlo move using Metropolis algorithm '''
    for i in range(N):
        for j in range(N):
                a = np.random.randint(N)
                b = np.random.randint(N)
                S0 =  lattice[a, b]
                #Sn = lattice[IN[a]%N,b] + lattice[a,IN[b]%N] + lattice[IP[a]%N,b] + lattice[a,IP[b]%N]
                Sn=lattice[(a - 1) % N, b] + lattice[(a+ 1) % N, b] + lattice[a, (b- 1) % N] + lattice[a, (b + 1) % N]
                dE = 2*S0*Sn
                if dE < 0 or np.random.random() < np.exp(-dE*beta):
                    lattice[i, j] = -lattice[i, j]  
                               
    return lattice
def calcEnergy(config,N):
    '''Energy of a given configuration'''
    energy = 0
    for i in range(len(config)):
        for j in range(len(config)):
            S = config[i,j]
            nb = config[(i+1)%N, j] + config[i,(j+1)%N] + config[(i-1)%N, j] + config[i,(j-1)%N]
            energy += -nb*S
    return energy/4.
def calcMag(lattice):
    '''Magnetization '''
    mag = np.sum(lattice)
    return mag  
for l in L:     
    n1, n2 = 1.0/(mcSteps), 1.0/(mcSteps*mcSteps)   
    E,M,C,X,BC = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)      
    for tt in range(nt):
        E1 = M1 = E2 = M2 = M3 = 0
        lattice = init_lattice(l)
        iT=1.0/T[tt]; iT2=iT*iT;
        
        for i in range(eqSteps):         # equilibrate
            mcmove(lattice, iT,l)           # Monte Carlo moves
        for i in range(mcSteps):
            mcmove(lattice, iT,l)  
            Ene = calcEnergy(lattice,l)   
            Mag = calcMag(lattice)        
            E1 = E1 + Ene
            M1 = M1 + Mag
            M2 = M2 + Mag*Mag 
            M3 = M3+Mag*Mag*Mag*Mag
            E2 = E2 + Ene*Ene               
    
        E[tt] = n1*E1
        M[tt] = n1*M1
        C[tt] = (n1*E2 - n2*E1*E1)*iT2
        X[tt] = (n1*M2 - n2*M1*M1)*iT
        BC[tt]= 1- (M3*n1)/ (3*M2*n1)
        if(T[tt]==1.53 or T[tt]==2.2541379310344825  or T[tt]==3.28):
            configurations.append(lattice)
    
    El.append (E)
    Ml.append (M)
    Cl.append (C)
    Xl.append (X)
    BCl.append (BC)
    print(l)
print(configurations)  

plt.xlabel("Temperature (T)")
plt.ylabel("Energy")
plt.plot(T,El[0], 'ro-',label='L=8', linewidth=0.5)
plt.plot(T,El[1], 'bs-',label='L=16', linewidth=0.5)
plt.plot(T,El[2],'g^-', label='L=35', linewidth=0.5)
plt.legend()
plt.savefig('../Ising/Figures/Energy.eps', format='eps')
plt.show()
plt.xlabel("Temperature (T)")
plt.ylabel("Magnetization ")
plt.plot(T,abs(Ml[0]), 'ro-',label='L=8', linewidth=0.5)
plt.plot(T,abs(Ml[1]), 'bs-',label='L=16', linewidth=0.5)
plt.plot(T,abs(Ml[2]),'g^-', label='L=35', linewidth=0.5)
plt.legend()
plt.savefig('../Ising/Figures/Magnetization.eps', format='eps')
plt.show()
plt.xlabel("Temperature (T)")
plt.ylabel("Specific Heat ")
plt.plot(T,Cl[0], 'ro-',label='L=8', linewidth=0.5)
plt.plot(T,Cl[1], 'bs-',label='L=16', linewidth=0.5)
plt.plot(T,Cl[2],'g^-', label='L=35', linewidth=0.5)
plt.legend()
plt.savefig('C:/Users/Gunel/Desktop/Ising/Figures/Specific_Heat.eps', format='eps')
plt.show()
plt.xlabel("Temperature (T)")
plt.ylabel("Susceptibility ")
plt.plot(T,Xl[0], 'ro-',label='L=8', linewidth=0.5)
plt.plot(T,Xl[1], 'bs-',label='L=16', linewidth=0.5)
plt.plot(T,Xl[2],'g^-', label='L=35', linewidth=0.5)
plt.legend()
plt.savefig('../Ising/Figures/Susceptibility.eps', format='eps')
plt.show()
plt.xlabel("Temperature (T)")
plt.ylabel("Binder’s cumulant ")
plt.plot(T,BCl[0], 'ro-',label='L=8', linewidth=0.5)
plt.plot(T,BCl[1], 'bs-',label='L=16', linewidth=0.5)
plt.plot(T,BCl[2],'g^-', label='L=35', linewidth=0.5)
plt.legend()
plt.savefig('../Ising/Figures/Binders_cumulant.eps', format='eps')
plt.show() 
plt.pcolormesh(configurations[0], cmap=plt.cm.YlGn);
plt.pcolormesh(configurations[1], cmap=plt.cm.YlGn);
plt.pcolormesh(configurations[2], cmap=plt.cm.YlGn);
plt.pcolormesh(configurations[3], cmap=plt.cm.YlGn);
plt.pcolormesh(configurations[4], cmap=plt.cm.YlGn);
plt.pcolormesh(configurations[5], cmap=plt.cm.YlGn);
plt.pcolormesh(configurations[6], cmap=plt.cm.YlGn);
plt.pcolormesh(configurations[7], cmap=plt.cm.YlGn);
plt.pcolormesh(configurations[8], cmap=plt.cm.YlGn);