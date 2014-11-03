#import libraries
import numpy as np
#import cmath
#import matplotlib as mpl
import matplotlib.pyplot as plt
import random

###############################start###########################################
def quantumwalk(steps=50,decotype='n',deco=0, returntype='p'): 
    #qubits |0> and |1>
    q0 = np.matrix('1;0')
    q1 = np.matrix('0;1')
    
    #variables changeable by user
    #steps = stepnumber #number of steps of the quantum walk to take
    #initial = q0 #initial coin is |0> - assume inital position is |0>
    initial = 1/np.sqrt(2)*np.matrix('1;1j') #balanced initial coin
    decoherence = deco #decoherence rate: probability of a decoherence event occuring per time step
    
    #set coin flip operator: I (tensor prod) H
    H = 1/np.sqrt(2)*np.matrix('1 1;1 -1') #hadamard matrix
    I_P = np.matrix(np.identity(2*steps+1))
    I_C = np.matrix(np.identity(2))
    coin_flip = np.kron(I_P,H)
    
    #set conditional translation operator, S
    S_0 = np.zeros((2*steps+1,2*steps+1))
    for x in range(0,2*steps):
        ket = np.zeros((2*steps+1,1))
        bra = np.zeros((2*steps+1,1))
        ket[x+1]=1
        bra[x]=1
        S_0 = S_0 + np.outer(ket,bra)
        
    S_1 = np.zeros((2*steps+1,2*steps+1))
    for x in range(1,2*steps+1):
        ket = np.zeros((2*steps+1,1))
        bra = np.zeros((2*steps+1,1))
        ket[x-1]=1
        bra[x]=1
        S_1 = S_1 + np.outer(ket,bra)
    
    S = np.kron(S_0,np.outer(q0,q0))+np.kron(S_1,np.outer(q1,q1))
    
    #set initial state of system
    system = np.zeros((4*steps+2,4*steps+2))    #empty density matrix of system: position density matrix (tensor product) coin density matrix
    system[2*steps:2*steps+2,2*steps:2*steps+2] = np.outer(np.conjugate(initial),initial)    #initial state of system (have to explicitly put conjugate into outer product)
    
    #Quantum Walk
    if decotype == 'n':
        for x in range(steps):
            system = S*coin_flip*system*coin_flip*np.matrix.getH(S) #regular quantum walk - no noise at all
    #Position decoherence version
    elif decotype == 'p':
        for x in range(steps):
            system_noiseless = S*coin_flip*system*coin_flip*np.matrix.getH(S) #step of the walk where there is no noise at that step
            system = (1-decoherence)*system_noiseless   #with probability (1-decoherence), the system is not affected by noise
            for i in range(2*steps+1):
                 ket_i = np.zeros((2*steps+1,1))
                 ket_i[i] = 1
                 M_i = np.outer(ket_i,ket_i)
                 M_P = np.kron(M_i,I_C)
                 system = system + decoherence*M_P*system_noiseless*np.matrix.getH(M_P) #with probability decoherence, the system's position is measured
    #Coin decoherence version
    elif decotype == 'c':
        for x in range(steps):
            system_noiseless = S*coin_flip*system*coin_flip*np.matrix.getH(S) #step of the walk where there is no noise at that step
            system = (1-decoherence)*system_noiseless   #with probability (1-decoherence), the system is not affected by noise
            for i in range(2):
                ket_i = np.zeros((2,1))
                ket_i[i] = 1
                M_i = np.outer(ket_i,ket_i)
                M_C = np.kron(I_P,M_i)
                system = system + decoherence*M_C*system_noiseless*np.matrix.getH(M_C) #with probability decoherence, the coin is measured
    #Coin and Position
    elif decotype == 'cp':
        for x in range(steps):
            system_noiseless = S*coin_flip*system*coin_flip*np.matrix.getH(S) #step of the walk where there is no noise at that step
            system = (1-decoherence)*system_noiseless   #with probability (1-decoherence), the system is not affected by noise
            for i in range(2*(2*steps+1)):
                ket_i = np.zeros((2*(2*steps+1),1))
                ket_i[i] = 1
                M_CP = np.outer(ket_i,ket_i)
                system = system + decoherence*M_CP*system_noiseless*np.matrix.getH(M_CP)
     #Quantum Walk
    elif decotype == 'h':
        for x in range(steps):
            ranvar = random.gauss(0.5,np.sqrt(deco)*0.01)
            H = 1/np.sqrt(2)*np.matrix([[np.sqrt(ranvar), np.sqrt(1-ranvar)], [np.sqrt(1-ranvar), -1*np.sqrt(ranvar)]])
            coin_flip = np.kron(I_P,H)
            system = S*coin_flip*system*coin_flip*np.matrix.getH(S) #regular quantum walk - no noise at all
                 
    else:
        print('Please pick one of the options')
       
    #Measurement
    p = np.zeros((2*steps+1,1))
    for i in range(2*steps+1):
        ket_i = np.zeros((2*steps+1,1))
        ket_i[i] = 1
        M_i = np.outer(ket_i,ket_i)
        M_P = np.kron(M_i,I_C)
        p[i]= np.trace(M_P*np.matrix.getH(M_P)*system)
    
    
    #Return type
    #plot
    if returntype == 'p':
        steplabels = range(-steps,steps+2,2)
        plt.plot(steplabels,p[0:2*steps+1:2])
    elif returntype == 's':
        return system
    elif returntype == 'm':
        return p
    elif returntype == 'v':
        devcopy = list(p[0:2*steps+1:2])
        steplabels2=list(range(-steps,steps+2,2))
        for i in range(0,51):
            devcopy[i]=devcopy[i]*steplabels2[i]
        return np.std(np.array(devcopy))
    elif returntype == 'a':
        return p
    else:
        steplabels = range(-steps,steps+2,2)
        plt.plot(steplabels,p[0:2*steps+1:2])      
    
    
    
       
#####################################end#######################################


