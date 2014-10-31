# -*- coding: utf-8 -*-
"""
Created on Sun Oct 26 19:44:22 2014

@author: Rosanna
"""
#import libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#qubits |0> and |1>
q0 = np.matrix('1;0')
q1 = np.matrix('0;1')

#variables changeable by user
steps = 50 #number of steps of the quantum walk to take
#initial = q0 #initial coin is |0> - assume inital position is |0>
initial = 1/np.sqrt(2)*np.matrix('1;1j') #balanced initial coin
decoherences = 15 #decoherence rate: probability of a decoherence event occuring per time step

#ANIMATION SET-UP
#set up figure
fig = plt.figure()
ax = plt.axes(xlim=(-steps, steps), ylim=(0, 0.12))
line, = ax.plot([], [], lw=2)
# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,
#set up x values (even numbers of steps)
steplabels = range(-steps,steps+2,2)

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

#THE QUANTUM WALK - WITHOUT NOISE
quantum_control = np.zeros((4*steps+2,4*steps+2))    #empty density matrix of system: position density matrix (tensor product) coin density matrix
quantum_control[2*steps:2*steps+2,2*steps:2*steps+2] = np.outer(np.conjugate(initial),initial)    #initial state of system (have to explicitly put conjugate into outer product)
for x in range(steps):
    quantum_control = S*coin_flip*quantum_control*coin_flip*np.matrix.getH(S) #regular quantum walk - no noise at all
#Measurement
p_quantum_control = np.zeros((2*steps+1,1))
for i in range(2*steps+1):
    ket_i = np.zeros((2*steps+1,1))
    ket_i[i] = 1
    M_i = np.outer(ket_i,ket_i)
    M_P = np.kron(M_i,I_C)
    p_quantum_control[i]= np.trace(M_P*np.matrix.getH(M_P)*quantum_control)
    
#THE QUANTUM WALK - WITH NOISE
p = np.zeros((2*steps+1,decoherences+1))
decoherences_list = np.zeros((decoherences+1,1))
stdevs = np.zeros((decoherences+1,1))
for j in range(decoherences+1):
    print(j)
    decoherence = np.cos((j)*np.pi*0.5/decoherences)
    decoherences_list[j] =  decoherence
    #set initial state of system
    system = np.zeros((4*steps+2,4*steps+2))    #empty density matrix of system: position density matrix (tensor product) coin density matrix
    system[2*steps:2*steps+2,2*steps:2*steps+2] = np.outer(np.conjugate(initial),initial)    #initial state of system (have to explicitly put conjugate into outer product)
    #The Quantum Walk
    for x in range(steps):
        system_noiseless = S*coin_flip*system*coin_flip*np.matrix.getH(S) #step of the walk where there is no noise at that step
        system = (1-decoherence)*system_noiseless   #with probability (1-decoherence), the system is not affected by noise
        for i in range(2*steps+1):
            ket_i = np.zeros((2*steps+1,1))
            ket_i[i] = 1
            M_i = np.outer(ket_i,ket_i)
            M_P = np.kron(M_i,I_C)
            system = system + decoherence*M_P*system_noiseless*np.matrix.getH(M_P) #with probability decoherence, the system's position is measured
    #Measurement
    for i in range(2*steps+1):
        ket_i = np.zeros((2*steps+1,1))
        ket_i[i] = 1
        M_i = np.outer(ket_i,ket_i)
        M_P = np.kron(M_i,I_C)
        p[i,j]= np.trace(M_P*np.matrix.getH(M_P)*system)
    stdevs[j] = np.std( p[0:2*steps+1:2,j])
    
def animate(j):
    y = p[0:2*steps+1:2,j]
    line.set_data(steplabels, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=decoherences+1, interval=80, blit=True)

plt.show()

#plot measurement probabilities
plt.figure()
plt.plot(decoherences_list,stdevs)