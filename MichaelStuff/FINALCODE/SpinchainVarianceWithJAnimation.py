#NEW QUANTUM WALK/SPIN CHAIN SCRIPT
##################################
#PREAMBLE
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
##################################
#INITIAL VARIABLES
STEPS = 10
POSITIONS = 2*STEPS+1
CHAINLENGTH = 3
results=[]
J = 0
J1 = 0.01
I_P = np.identity((POSITIONS))
I_C = np.identity((2))
##################################
#HADAMARD OPERATOR
H = (1/np.sqrt(2))*np.matrix([[1,  1],
                              [1, -1]])
H = np.kron(I_P,H)
for z in range(CHAINLENGTH-1):
    H = np.kron(H,I_C)
##################################
#SHIFT OPERATOR
SHIFTUP = np.matrix(np.zeros((POSITIONS,POSITIONS)))
for y in range(POSITIONS-1):
    ket = np.matrix(np.zeros((POSITIONS,1)))
    bra = np.matrix(np.zeros((POSITIONS,1)))
    ket[y+1]=1
    bra[y]=1
    SHIFTUP = SHIFTUP + np.outer(ket,bra)
COINUP = np.matrix([[1, 0],
                    [0, 0]])
SHIFTUP = np.kron(SHIFTUP,COINUP)
    
SHIFTDOWN = np.matrix(np.zeros((POSITIONS,POSITIONS)))
for x in range(1,POSITIONS):
    ket = np.matrix(np.zeros((POSITIONS,1)))
    bra = np.matrix(np.zeros((POSITIONS,1)))
    ket[x-1]=1
    bra[x]=1
    SHIFTDOWN = SHIFTDOWN + np.outer(ket,bra)
COINDOWN = np.matrix([[0, 0],
                      [0, 1]])
SHIFTDOWN = np.kron(SHIFTDOWN,COINDOWN)

for w in range(CHAINLENGTH-1):
    SHIFTUP = np.kron(SHIFTUP,I_C)
    SHIFTDOWN = np.kron(SHIFTDOWN,I_C)
SHIFT = SHIFTUP+SHIFTDOWN

##################################
#SPINCHAIN OPERATOR
pauli_x = np.matrix(('0 1;1 0'))
pauli_y = np.matrix(('0 -1j;1j 0'))
pauli_z = np.matrix(('1 0;0 -1'))
pauli_i = np.matrix(('1 0;0 1'))

CHAIN1 = 0
CHAIN2 = 0
for v in range(CHAINLENGTH-1):
    COMPONENT1 = 1
    COMPONENT2 = 1
    for u in range(CHAINLENGTH):
        if u==v:
            COMPONENT1 = np.kron(COMPONENT1,pauli_x)
            COMPONENT2 = np.kron(COMPONENT2,pauli_y)
        elif u==v+1:
            COMPONENT1 = np.kron(COMPONENT1,pauli_x)
            COMPONENT2 = np.kron(COMPONENT2,pauli_y)
        else:
            COMPONENT1 = np.kron(COMPONENT1,pauli_i)
            COMPONENT2 = np.kron(COMPONENT2,pauli_i)

    if v==0:
        CHAIN1 = COMPONENT1 + COMPONENT2
    else:
        CHAIN2 = CHAIN2 + COMPONENT1 + COMPONENT2

#####################################
#CLEAR MEMORY
del SHIFTUP
del SHIFTDOWN
del COINDOWN
del COINUP
del z
del y
del ket
del bra
del x
del w
del pauli_x
del pauli_y
del pauli_z
del pauli_i
del v
del u
del COMPONENT1
del COMPONENT2
######################################
#INITIAL STATE
SYSTEM = np.matrix(np.zeros((POSITIONS,1)))
SYSTEM[STEPS] = (1)
SYSTEM = np.outer(SYSTEM,SYSTEM)
for ab in range(CHAINLENGTH):
    SYSTEM = np.kron(SYSTEM,I_C)
    
INITIALCOIN =               np.matrix([[1],
                                       [0]])
INITIALCOIN = np.outer(INITIALCOIN,INITIALCOIN)
INITIALCOIN = np.kron(I_P,INITIALCOIN)
for aC in range(CHAINLENGTH-1):
    INITIALCOIN = np.kron(INITIALCOIN,I_C)    

INITIALCHAIN = (1/np.sqrt(2))*np.matrix([[1],
                          [1]])
INITIALCHAIN1 = np.outer(INITIALCHAIN,INITIALCHAIN)
INITIALCHAIN = np.kron(I_P,I_C)        
for ad in range(CHAINLENGTH-1):
    INITIALCHAIN = np.kron(INITIALCHAIN,INITIALCHAIN1)   

SYSTEM = np.matrix(SYSTEM)*np.matrix(INITIALCOIN)*np.matrix(INITIALCHAIN)
SYSTEM0 = SYSTEM

del INITIALCHAIN
del INITIALCOIN
del INITIALCHAIN1
##################################

for decoherences in range(0,100):
    J=decoherences/100
    SYSTEM = SYSTEM0
############################################################ 
#LAST BIT OF CHAIN CODE
    CHAIN = J1*CHAIN1 + J*CHAIN2
    EXP1 = linalg.expm2(-1j*CHAIN)
    EXP2 = linalg.expm2(1j*CHAIN)
    EXP1 = np.kron(I_P,EXP1)
    EXP2 = np.kron(I_P,EXP2)
###############################################################
    #QUANTUM WALK
    for r in range(STEPS):
        SYSTEM = np.matrix(EXP1)*np.matrix(SYSTEM)*np.matrix(EXP2)        
        SYSTEM = np.matrix(H)*np.matrix(SYSTEM)*np.matrix(np.matrix.getH(H))
        SYSTEM = np.matrix(SHIFT)*np.matrix(SYSTEM)*np.matrix(np.matrix.getH(SHIFT))
    ##################################
    #MEASUREMENT
    p = np.zeros((POSITIONS,1),dtype=complex)
    for i in range(POSITIONS):
        ket_i = np.zeros((POSITIONS,1))
        ket_i[i] = 1
        M_i = np.outer(ket_i,ket_i)
        M_P = np.kron(M_i,I_C)
        for q in range(CHAINLENGTH-1):
            M_P = np.kron(M_P,I_C)
        p[i]= np.trace(M_P*SYSTEM*np.matrix.getH(M_P))
        
    COPY = list(np.absolute(p[0:POSITIONS:2]))
    LABELS=list(range(-STEPS,STEPS+2,2))
    ##################################
    #ANALYSIS
    variance = 0
    mean = 0
    for i in range(0,STEPS+1):
        mean = mean + COPY[i]*LABELS[i]
    mean = mean/(STEPS+1)
    for j in range(0,STEPS+1):
        variance = variance + COPY[j]*(LABELS[j] - mean)*(LABELS[j] - mean)
    results.append(variance)
    ##################################
plt.plot(results)
    
