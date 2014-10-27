# -*- coding: utf-8 -*-
"""
Created on Sun Oct 26 19:44:22 2014

@author: Rosanna
"""
#import libraries
import numpy as np
#import cmath
#import matplotlib as mpl
import matplotlib.pyplot as plt

#qubits |0> and |1>
q0 = np.matrix('1;0')
q1 = np.matrix('0;1')

#variables changeable by user
steps = 1000 #number of steps of the quantum walk to take
#initial = q0 #initial coin is |0> - assume inital position is |0>
initial = 1/np.sqrt(2)*np.matrix('1;1j') #balanced initial coin

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

#THE QUANTUM WALK
for x in range(steps):
    system = S*coin_flip*system*coin_flip*np.matrix.getH(S)

#pcolor(S)  #plots system matrix on grid so we can see non-zero values

#Measurement
p = np.zeros((2*steps+1,1))
for i in range(2*steps+1):
    ket_i = np.zeros((2*steps+1,1))
    ket_i[i] = 1
    M_i = np.outer(ket_i,ket_i)
    M_P = np.kron(M_i,I_C)
    p[i]= np.trace(M_P*np.matrix.getH(M_P)*system)
    
#plot measurement probabilities
steplabels = range(-steps,steps+2,2)
plt.plot(steplabels,p[0:2*steps+1:2])