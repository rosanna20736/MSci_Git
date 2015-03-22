#import libraries
import numpy as np
#import cmath
#import matplotlib as mpl
import matplotlib.pyplot as plt
import random

############################### random wieghted choice function ###########################################
def weighted_choice(weights):
    random.seed()
    totals = []
    running_total = 0

    for w in weights:
        running_total += w
        totals.append(running_total)

    rnd = random.random() * running_total
    for i, total in enumerate(totals):
        if rnd < total:
            return i

############################### variance measuring function ###########################################
def variance(p,steps):
    if np.count_nonzero(p) % 2:
        devcopy = list(p[0:2*steps+1:2])
        steplabels2=list(range(-steps,steps+2,2))
        varno = steps+1
    else:
        devcopy = list(p[1:2*steps+1:2])
        steplabels2=list(range(-steps+1,steps+1,2))
        varno = steps
    mean = 0
    variance = 0
    for i in range(0,varno):
        mean = mean + devcopy[i]*steplabels2[i]
    mean = mean/(varno)
    for j in range(0,varno):
        variance = variance + devcopy[j]*(steplabels2[j] - mean)*(steplabels2[j] - mean)
    return variance
    
################### measurement function ###############################################
def measure(system, steps):
    p = np.zeros((2*steps+1,1))
    I_C = np.matrix(np.identity(2))
    for i in range(2*steps+1):
        ket_i = np.zeros((2*steps+1,1))
        ket_i[i] = 1
        M_i = np.outer(ket_i,ket_i)
        M_P = np.kron(M_i,I_C)
        p[i]= np.trace(np.matrix.getH(M_P)*M_P*system)
    return p

############################### QUANTUM WALK FUNCTION ###########################################
#Function 'quantumwalk(steps, decoherence type, decoherence value, return type)
#Steps: 50 (Default)
#Decoherence Type: 'n' (none)(Default), 'c' (coin), 'p' (position), 'cp' (coin and position), 'H' (Hadamard), 'm' (Type of measurement chosen at random)
#Decoherence Value: 0 (Default)
#Return Type: 'p' (plot)(Default), 's' (system), 'm' (measurement), 'v' (variance), 'vb' (variance behaviour)
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
    if returntype == 'vb':
        variances=[]
    
    if decotype == 'n':
        for x in range(steps):
            system = S*coin_flip*system*coin_flip*np.matrix.getH(S) #regular quantum walk - no noise at all
            if returntype == 'vb':
                variances.append(variance(measure(system,steps),steps))
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
            if returntype == 'vb':
                variances.append(variance(measure(system,steps),steps))
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
            if returntype == 'vb':
                variances.append(variance(measure(system,steps),steps))
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
            if returntype == 'vb':
                variances.append(variance(measure(system,steps),steps))
    #Type of measurement randomly chosen if the system is measured
    elif decotype == 'm':
        p = 0.6
        for x in range(steps):
            noise_type = weighted_choice([p,(1-p)/(1+p),p*(1-p)/(1+p)]) #chooses position, coin, or both with weights given by p for position measurement assuming c and p indep.
            system_noiseless = S*coin_flip*system*coin_flip*np.matrix.getH(S) #step of the walk where there is no noise at that step
            system = (1-decoherence)*system_noiseless   #with probability (1-decoherence), the system is not affected by noise
            if noise_type == 0: #position measured
                for i in range(2*steps+1):
                    ket_i = np.zeros((2*steps+1,1))
                    ket_i[i] = 1
                    M_i = np.outer(ket_i,ket_i)
                    M_P = np.kron(M_i,I_C)
                    system = system + decoherence*M_P*system_noiseless*np.matrix.getH(M_P) #with probability decoherence, the system's position is measured 
            elif noise_type == 1: #coin measured
                for i in range(2):
                    ket_i = np.zeros((2,1))
                    ket_i[i] = 1
                    M_i = np.outer(ket_i,ket_i)
                    M_C = np.kron(I_P,M_i)
                    system = system + decoherence*M_C*system_noiseless*np.matrix.getH(M_C) #with probability decoherence, the coin is measured
            elif noise_type == 2: #both measured
                for i in range(2*(2*steps+1)):
                    ket_i = np.zeros((2*(2*steps+1),1))
                    ket_i[i] = 1
                    M_CP = np.outer(ket_i,ket_i)
                    system = system + decoherence*M_CP*system_noiseless*np.matrix.getH(M_CP)
            else:
                    print('none of coin/position/both chosen')
            if returntype == 'vb':
                variances.append(variance(measure(system,steps),steps))
    #Imperfect Hadamard Quantum Walk
    elif decotype == 'H':
        for x in range(steps):
            random.seed()
            noisy_theta = random.gauss(0.5*np.pi,np.sqrt(deco)*0.25*np.pi)
            noisy_H = np.matrix([[np.cos(noisy_theta/2),np.sin(noisy_theta/2)],[np.cos(noisy_theta/2),-np.sin(noisy_theta/2)]])
            noisy_coin_flip = np.kron(I_P,noisy_H)
            system = S*noisy_coin_flip*system*np.matrix.getH(noisy_coin_flip)*np.matrix.getH(S)
            if returntype == 'vb':
                p = measure(system,steps)
                variances.append(variance(p,steps))
    elif decotype == 'f':
        for x in range(steps):
            variances=(1-1/np.sqrt(2))*steps**2*(1-np.sqrt(2)/6*decoherence*steps+decoherence*(np.sqrt(2)-1))
        return variances
    else:
        print('Please pick one of the options')
       
    #Measurement
    p = measure(system,steps)
    
    #Return type
    #plot
    if returntype == 'p':
        steplabels = range(-steps,steps+2,2)
        plt.plot(steplabels,p[0:2*steps+1:2],label='p ='+ str(round(decoherence,2)))
    elif returntype == 's':
        return system
    elif returntype == 'm':
        return p
    elif returntype == 'v':
        return variance(p,steps)
    elif returntype == 'vb':
        varfit = np.polyfit(np.log(np.arange(1,steps+1)),np.log(variances),1,full=False)
        #print(varfit)
        return varfit[0]
    else:
        steplabels = range(-steps,steps+2,2)
        plt.plot(steplabels,p[0:2*steps+1:2])      
    
    
    
       
#####################################end#######################################


