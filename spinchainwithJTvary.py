#import libraries
import numpy as np
import matplotlib.pyplot as plt
import spinchaingenerator

#qubits |0> and |1>
q0 = np.matrix('1;0')
q1 = np.matrix('0;1')

chainlength = 2

#variables changeable by user
steps = 20 #number of steps of the quantum walk to take
JT = 2
#initial = q0 #initial coin is |0> - assume inital position is |0>
initial = 1/np.sqrt(2)*np.matrix('1;1j') #balanced initial coin
initialchain = np.matrix('0;1')
##set up x values (even numbers of steps)
steplabels = range(-steps,steps+2,2)
results = list(range(0,1000))

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


#set up intial system
initial1 = np.outer(np.conjugate(initial),initial)
system_AD = np.zeros((4*steps+2,4*steps+2))
system_AD[2*steps:2*steps+2,2*steps:2*steps+2] = initial1
initialchain = np.outer(np.conjugate(initialchain),initialchain)
for f in range(chainlength-1):
    system_AD = np.kron(system_AD,initialchain)
    
for x in range(0,1000):
    JT = x/10000
#    print(JT)
    p_AD = np.zeros(2*steps+1)
    #get spinchain
    spinchain = spinchaingenerator.U_spinchain('xy',JT)
    spinchain = np.kron(I_P,spinchain)
    
    #make overall operator
    OPERATOR1 = S*coin_flip*spinchain
    OPERATOR2 = np.matrix.getH(spinchain)*coin_flip*np.matrix.getH(S)
    
    #The Quantum Walk
    
    for y in range(steps):
        system_AD = OPERATOR1*system_AD*OPERATOR2
    
    
    
    #Measurement
    
    for i in range(0,2*steps+1,2):
        ket_i = np.zeros((2*steps+1,1))
        ket_i[i] = 1
        M_i = np.outer(ket_i,ket_i)
        M_P = np.kron(M_i,I_C)
        for z in range(chainlength-1):
            M_P = np.kron(M_P,I_C)
        p_AD[i]= np.trace(M_P*np.matrix.getH(M_P)*system_AD)
    
            
#    devcopy = list(p_AD[0:2*steps+1:2])
#    steplabels2=list(range(-steps,steps+2,2))
#    for i in range(0,steps+1):
#       devcopy[i]=devcopy[i]*steplabels2[i]
    mean = 0
    devcopy = 0
    steplabels2 = 0
    variance = 0
    devcopy = list(p_AD[0:2*steps+1:2])
    steplabels2=list(range(-steps,steps+2,2))
    for i in range(0,21):
        mean = mean + devcopy[i]*steplabels2[i]
    for j in range(0,21):
        variance = variance + devcopy[j]*(steplabels2[j] - mean)*(steplabels2[j] - mean)
    results[x] = variance
    print(variance)
    print(JT)
    


plt.plot(results)
plt.show()