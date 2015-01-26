#import libraries
import numpy as np
import matplotlib.pyplot as plt
import spinchain as spinchaingenerator

#qubits |0> and |1>
q0 = np.matrix('1;0')
q1 = np.matrix('0;1')

#variables changeable by user
chainlength = 5 #number of particles in spin chain
steps = 40 #number of steps of the quantum walk to take
#initial = q0 #initial coin is |0> - assume inital position is |0>
initial = 1/np.sqrt(2)*np.matrix('1;1j') #balanced initial coin
initialchain = q0

##set up x values (even numbers of steps)
steplabels = range(-steps,steps+2,2)

#set coin flip operator: I (tensor prod) H
H = 1/np.sqrt(2)*np.matrix('1 1;1 -1') #hadamard matrix
I_P = np.matrix(np.identity(2*steps+1))
I_C = np.matrix(np.identity(2))
coin_flip = np.kron(I_P,H)
for c in range(chainlength-1):
    coin_flip = np.kron(coin_flip,I_C)

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

S_0 = np.kron(S_0,np.outer(q0,q0))
S_1 = np.kron(S_1,np.outer(q1,q1))
for d in range(chainlength-1):
    S_0 = np.kron(S_0,I_C)
    S_1 = np.kron(S_1,I_C)

S = S_0 + S_1

#set up measurement array    
p = np.zeros(2*steps+1)

#set up intial system
system = np.zeros((4*steps+2,4*steps+2))
system[2*steps:2*steps+2,2*steps:2*steps+2] = np.outer(np.conjugate(initial),initial)
initialchain = np.outer(np.conjugate(initialchain),initialchain)
for f in range(chainlength-1):
    system = np.kron(system,initialchain)
    
#get spinchain
spinchain = spinchaingenerator.U_spinchain('xy')
spinchain = np.kron(I_P,spinchain)

#make overall operator - operator2 is hermitian conjugate of operator1
OPERATOR1 = S*coin_flip*spinchain
OPERATOR2 = np.matrix.getH(spinchain)*coin_flip*np.matrix.getH(S)

#The Quantum Walk
TIME = 0
for x in range(steps):
    system = OPERATOR1*system*OPERATOR2
    print(TIME)
    TIME = TIME+1

#Measurement
TIME2 = 0
for i in range(2*steps+1):
    ket_i = np.zeros((2*steps+1,1))
    ket_i[i] = 1
    M_i = np.outer(ket_i,ket_i)
    M_P = np.kron(M_i,I_C)
    for z in range(chainlength-1):
        M_P = np.kron(M_P,I_C)
    p[i]= np.trace(M_P*np.matrix.getH(M_P)*system)
    print(TIME2)
    TIME2 = TIME2+1

plt.plot(steplabels, p[0:2*steps+1:2])
plt.show()
