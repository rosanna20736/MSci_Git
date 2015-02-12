#SIMULATES QUANTUM WALK FOR AMPLITUDE AND PHASE DAMPING AND PLOTS VARIANCE WITH DECOHERENCE PROBABILITY

#import libraries
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as animation

plt.close('all')
#qubits |0> and |1>
q0 = np.matrix('1;0')
q1 = np.matrix('0;1')

#variables changeable by user
steps = 100 #number of steps of the quantum walk to take
#initial = q0 #initial coin is |0> - assume inital position is |0>
initial = 1/np.sqrt(2)*np.matrix('1;1j') #balanced initial coin
decoherences = 31 #decoherence rate: probability of a decoherence event occuring per time step

##ANIMATION SET-UP
##set up figure
#fig = plt.figure()
#ax = plt.axes(xlim=(-steps, steps), ylim=(0, 0.10))
#line, = ax.plot([], [], lw=2)
#line_AD, = ax.plot([], [], lw=2)
## initialization function: plot the background of each frame
#def init():
#   line.set_data([], [])
#   return line,
   
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
    
#THE QUANTUM WALK - WITH NOISE
p_AD = np.zeros((2*steps+1,decoherences+1))
p_PD = np.zeros((2*steps+1,decoherences+1))
decoherences_list = []
variances_AD = []
variances_PD = []
for j in range(decoherences+1):
    print(j)
    decoherence = 1 - np.log(1+(np.e-1)*j/decoherences)
    decoherences_list.append(decoherence)
    #set initial state of system
    system_AD = np.zeros((4*steps+2,4*steps+2))
    system_AD[2*steps:2*steps+2,2*steps:2*steps+2] = np.outer(np.conjugate(initial),initial)
    system_PD = np.zeros((4*steps+2,4*steps+2))
    system_PD[2*steps:2*steps+2,2*steps:2*steps+2] = np.outer(np.conjugate(initial),initial)
    #The Quantum Walk
    for x in range(steps):
        constant=decoherence
        E0 = np.matrix([[1, 0], [0, np.sqrt(1-constant)]]) #E1 is the same for both amplitude and phase damping
        E1 = np.matrix([[0, np.sqrt(constant)], [0, 0]])
        E1_PD = np.matrix([[0, 0], [0, np.sqrt(constant)]])
        E0 = np.kron(I_P,E0)
        E1 = np.kron(I_P,E1)
        E1_PD = np.kron(I_P,E1_PD)
        system_AD = S*coin_flip*E0*system_AD*np.matrix.getH(E0)*coin_flip*np.matrix.getH(S) + S*coin_flip*E1*system_AD*np.matrix.getH(E1)*coin_flip*np.matrix.getH(S)
        system_PD = S*coin_flip*E0*system_PD*np.matrix.getH(E0)*coin_flip*np.matrix.getH(S) + S*coin_flip*E1_PD*system_PD*np.matrix.getH(E1_PD)*coin_flip*np.matrix.getH(S)
        
    #Measurement
    for i in range(2*steps+1):
        ket_i = np.zeros((2*steps+1,1))
        ket_i[i] = 1
        M_i = np.outer(ket_i,ket_i)
        M_P = np.kron(M_i,I_C)
        p_AD[i,j]= np.trace(M_P*np.matrix.getH(M_P)*system_AD)
        p_PD[i,j]= np.trace(M_P*np.matrix.getH(M_P)*system_PD)
#    plt.figure(1)
#    plt.plot(steplabels,p_AD[0:2*steps+1:2,j],label='$\gamma$ ='+ str(round(decoherence,2)))
#    plt.figure(2)
#    plt.plot(steplabels,p_PD[0:2*steps+1:2,j],label='$\gamma$ ='+ str(round(decoherence,2)))

    devcopy_AD = list(p_AD[0:2*steps+1:2,j])
    devcopy_PD = list(p_PD[0:2*steps+1:2,j])
    steplabels2=list(range(-steps,steps+2,2))
    mean_AD = 0
    variance_AD = 0
    mean_PD = 0
    variance_PD = 0
    for i in range(0,steps+1):
        mean_AD = mean_AD + devcopy_AD[i]*steplabels2[i]
        mean_PD = mean_PD + devcopy_PD[i]*steplabels2[i]
    mean_AD = mean_AD/(steps+1)
    mean_PD = mean_PD/(steps+1)
    for j in range(0,steps+1):
        variance_AD = variance_AD + devcopy_AD[j]*(steplabels2[j] - mean_AD)*(steplabels2[j] - mean_AD)
        variance_PD = variance_PD + devcopy_PD[j]*(steplabels2[j] - mean_PD)*(steplabels2[j] - mean_PD)
    variances_AD.append(variance_AD)
    variances_PD.append(variance_PD)
    
#plt.figure(1)
#plt.xlabel('particle position, x',fontsize='large')
#plt.ylabel('probability particle is found at position x',fontsize='large')
#plt.axis([-steps,steps,0,np.max(np.max(p_AD))])
#plt.legend()
#plt.show
#
#plt.figure(2)
#plt.xlabel('particle position, x',fontsize='large')
#plt.ylabel('probability particle is found at position x',fontsize='large')
#plt.axis([-steps,steps,0,np.max(np.max(p_PD))])
#plt.legend()
#plt.show

plt.figure(3)
plt.plot(decoherences_list,variances_AD,label='amplitude damping')
plt.plot(decoherences_list,variances_PD,label='phase damping')
plt.xlabel('probability of decoherence event',fontsize='large')
plt.ylabel('variance of walk',fontsize='large')
plt.legend()
plt.axis([0,1,0,np.max([variances_AD,variances_PD])])
plt.show

#def animate(j): 
##    line.set_data(steplabels, p[0:2*steps+1:2,j])
##    line_C.set_data(steplabels, p_C[0:2*steps+1:2,j])
##    line_CP.set_data(steplabels, p_CP[0:2*steps+1:2,j])
#    line_AD.set_data(steplabels, p_AD[0:2*steps+1:2,j])
#    #return line, line_C, line_CP, line_AD,
#    return line_AD,
## call the animator.  blit=True means only re-draw the parts that have changed.
#anim = animation.FuncAnimation(fig, animate, init_func=init, frames=decoherences+1, interval=80, blit=True)
#plt.show()
