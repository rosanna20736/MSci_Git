#NEW QUANTUM WALK/SPIN CHAIN SCRIPT
##################################
#PREAMBLE
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
##################################
#INITIAL VARIABLES
STEPS = 150
POSITIONS = 2*STEPS+1
CHAINLENGTH = 2
results=[]
J = 0.1
J1 =0.
I_P = np.identity((POSITIONS))
I_C = np.identity((2))
q0 = np.matrix('1;0')
q1 = np.matrix('0;1')
    
###################################
#HADAMARD OPERATOR
H = (1/np.sqrt(2))*np.matrix([[1,  1],
                              [1, -1]])
H = np.kron(I_P,H)
for z in range(CHAINLENGTH-1):
    H = np.kron(H,I_C)
#print(np.matrix(H)*np.matrix.getH(np.matrix(H)))
##################################
#SHIFT OPERATOR
#set conditional translation operator, S
S_0 = np.zeros((4*STEPS+2,4*STEPS+2))
for x in range(0,2*STEPS):
    ket = np.zeros((2*STEPS+1,1))
    bra = np.zeros((2*STEPS+1,1))
    ket[x+1]=1
    bra[x]=1
    S_TEMP = np.kron(np.outer(ket,bra),np.outer(q0,q0))
    S_0 = S_0 + S_TEMP
#    S_0 = S_0 + np.outer(ket,bra)
    
S_1 = np.zeros((4*STEPS+2,4*STEPS+2))
for x in range(1,2*STEPS+1):
    ket = np.zeros((2*STEPS+1,1))
    bra = np.zeros((2*STEPS+1,1))
    ket[x-1]=1
    bra[x]=1
    S_TEMP1 = np.kron(np.outer(ket,bra),np.outer(q1,q1))
    S_1 = S_1 + S_TEMP1
#    S_1 = S_1 + np.outer(ket,bra)

#S = np.kron(S_0,np.outer(q0,q0))+np.kron(S_1,np.outer(q1,q1))
S = S_1+S_0
#print(np.matrix(S)*np.matrix.getH(np.matrix(S)))
S = np.kron(S,I_C)

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
A = 0.5*(J*CHAIN1 + J*CHAIN2)
A = 1j*A
EXP1 = linalg.expm(-1*A)
EXP2 = linalg.expm(A)
EXP1 = np.kron(I_P,EXP1)
EXP2 = np.kron(I_P,EXP2)
######################################
#INITIAL STATE   
INITIALCOIN =               (1/np.sqrt(2))*np.matrix([[1],
                                       [1]])
INITIALCOIN = np.outer(INITIALCOIN,INITIALCOIN)
#INITIALCHAIN = np.matrix([[0],
#                         [1]])
#INITIALCHAIN = np.outer(INITIALCHAIN,INITIALCHAIN)
INITIALCHAIN = np.matrix([[0,  0],
                                         [0,  1]])  
#SYSTEM = np.zeros((4*STEPS+2,4*STEPS+2))
#SYSTEM[2*STEPS:2*STEPS+2,2*STEPS:2*STEPS+2] = np.outer(np.conjugate(INITIALCOIN),INITIALCOIN)
#for ad in range(CHAINLENGTH-1):
#    SYSTEM = np.kron(SYSTEM,INITIALCHAIN)
SYSTEM = np.zeros((2*STEPS+1))
SYSTEM[STEPS+1] = 1
SYSTEM = np.outer(SYSTEM,SYSTEM)


A = np.kron(INITIALCOIN,INITIALCHAIN)
SYSTEM = np.kron(SYSTEM,A)

######################################
#INITIAL STATE
INITIALCOIN1 =               (1/np.sqrt(2))*np.matrix([[1],
                                       [-1]])
INITIALCOIN1 = np.outer(INITIALCOIN1,INITIALCOIN1)
#INITIALCHAIN = np.matrix([[0],
#                         [1]])
#INITIALCHAIN = np.outer(INITIALCHAIN,INITIALCHAIN)
INITIALCHAIN1 = np.matrix([[0,  0],
                                         [0,  1]])  
#SYSTEM = np.zeros((4*STEPS+2,4*STEPS+2))
#SYSTEM[2*STEPS:2*STEPS+2,2*STEPS:2*STEPS+2] = np.outer(np.conjugate(INITIALCOIN),INITIALCOIN)
#for ad in range(CHAINLENGTH-1):
#    SYSTEM = np.kron(SYSTEM,INITIALCHAIN)
SYSTEM1 = np.zeros((2*STEPS+1))
SYSTEM1[STEPS+1] = 1
SYSTEM1 = np.outer(SYSTEM1,SYSTEM1)


A1 = np.kron(INITIALCOIN1,INITIALCHAIN1)
SYSTEM1 = np.kron(SYSTEM1,A1)

########################################
#QUANTUM WALK
p = np.zeros((POSITIONS,STEPS),dtype=complex)
for r in range(STEPS):
  
#    SYSTEM = np.matrix(H)*np.matrix(SYSTEM)*np.matrix(np.matrix.getH(H))
#    SYSTEM = np.matrix(S)*np.matrix(SYSTEM)*np.matrix(np.matrix.getH(S))
#    SYSTEM = np.matrix(EXP2)*np.matrix(SYSTEM)*np.matrix(EXP1)     
##
#    SYSTEM2 = np.matrix(H)*np.matrix(SYSTEM1)*np.matrix(np.matrix.getH(H))
#    SYSTEM2 = np.matrix(S)*np.matrix(SYSTEM1)*np.matrix(np.matrix.getH(S))
#    SYSTEM2 = np.matrix(EXP1)*np.matrix(SYSTEM1)*np.matrix(EXP2)        
#    system = SYSTEM - SYSTEM2
#    system = np.matrix.conj(system)*np.matrix(system)
##################################
    #MEASUREMENT
#    x = 0
#    for i in range(POSITIONS):
#        ket_i = np.zeros((POSITIONS,1))
#        ket_i[i] = 1
#        M_i = np.outer(ket_i,ket_i)
#        M_P = np.kron(M_i,I_C)
#        for q in range(CHAINLENGTH-1):
#            M_P = np.kron(M_P,I_C)
#        p[i,r]= (np.trace(M_P*system*np.matrix.getH(M_P)))
#    for j in range(POSITIONS):
#        x = x + np.absolute(p[j,r])
#    results.append(0.5*x)
    
    #Coin Meaurement
    x = np.zeros((2,2))
    y = np.zeros((2,2))
    for i in range(POSITIONS):
        ket_i = np.zeros((POSITIONS,1))
        ket_i[i] = 1
        M_i = np.matrix.getH(ket_i)
        M_P = np.kron(M_i,I_C)
        for j in range(2):
            ket_j = np.zeros(2)
            ket_j[j] = 1
            M_z = np.matrix.getH(ket_j)
            M_C = np.kron(M_P,M_z) 
        x = x + np.matrix(M_C)*SYSTEM*np.matrix.getH(M_C)
        y = y + np.matrix(M_C)*SYSTEM1*np.matrix.getH(M_C)
    z = linalg.sqrtm(np.conj(x-y)*(x-y))
#    z = np.conj(x-y)*(x-y)
    results.append(0.5*np.absolute(np.trace(z)))
#    results.append(0.5*(np.absolute(np.trace(MESUP*system*np.matrix.getH(MESUP)))+np.absolute(np.trace(MESDOWN*system*np.matrix.getH(MESDOWN)))))
    
    SYSTEM = np.matrix(H)*np.matrix(SYSTEM)*np.matrix(np.matrix.getH(H))
    SYSTEM = np.matrix(S)*np.matrix(SYSTEM)*np.matrix(np.matrix.getH(S))
    SYSTEM = np.matrix(EXP2)*np.matrix(SYSTEM)*np.matrix(EXP1)     
#
    SYSTEM1 = np.matrix(H)*np.matrix(SYSTEM1)*np.matrix(np.matrix.getH(H))
    SYSTEM1 = np.matrix(S)*np.matrix(SYSTEM1)*np.matrix(np.matrix.getH(S))
    SYSTEM1 = np.matrix(EXP1)*np.matrix(SYSTEM1)*np.matrix(EXP2)           
    
     
plt.plot(results)

