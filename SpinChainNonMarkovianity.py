#NEW QUANTUM WALK/SPIN CHAIN SCRIPT
##################################
#PREAMBLE
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
##################################
#INITIAL VARIABLES
STEPS = 23
POSITIONS = 2*STEPS+1
CHAINLENGTH =5
results=[]
results2=[]
J = 0.2
J1 =0.1
I_P = np.identity((POSITIONS))
I_C = np.identity((2))
q0 = np.matrix('1;0')
q1 = np.matrix('0;1')
    
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
A = (J1*CHAIN1 + J*CHAIN2)
A = 1j*A
EXP1 = linalg.expm(-1*A)
EXP2 = linalg.expm(A)
######################################
#INITIAL STATE   
#INITIALCOIN =               (1/np.sqrt(2))*np.matrix([[1],
#                                       [1]])
INITIALCOIN =               np.matrix([[1],
                                       [0]])
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
A = np.kron(INITIALCOIN,INITIALCHAIN)                                         
for x in range(CHAINLENGTH-2):
    A = np.kron(A,INITIALCHAIN)
SYSTEM = A
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
A1 = np.kron(INITIALCOIN1,INITIALCHAIN1)                                         
for x in range(CHAINLENGTH-2):
    A1 = np.kron(A1,INITIALCHAIN1)

SYSTEM1 = A1
########################################
#QUANTUM WALK

for r in range(STEPS):
    
    #Coin Meaurement
    x = np.zeros((2,2))
    y = np.zeros((2,2))
    for i in range(2):
        ket_i = np.zeros(2)
        ket_i[i] = 1
        M_i = np.matrix.getH(ket_i)
        M_A= np.kron(I_C,M_i)
        if CHAINLENGTH == 2:
            x = x + np.matrix(M_A)*SYSTEM*np.matrix.getH(M_A)
        if CHAINLENGTH > 2:
            for i1 in range(2):
                ket_i = np.zeros(2)
                ket_i[i1] = 1
                M_i = np.matrix.getH(ket_i)
                M_B = np.kron(M_A,M_i)
                if CHAINLENGTH == 3:
                    x = x + np.matrix(M_B)*SYSTEM*np.matrix.getH(M_B)
                if CHAINLENGTH > 3:
                    for i2 in range(2):
                        ket_i = np.zeros(2)
                        ket_i[i2] = 1
                        M_i = np.matrix.getH(ket_i)
                        M_C = np.kron(M_B,M_i)
                        if CHAINLENGTH == 4:
                            x = x + np.matrix(M_C)*SYSTEM*np.matrix.getH(M_C)
                        if CHAINLENGTH > 4:
                            for i3 in range(2):
                                ket_i = np.zeros(2)
                                ket_i[i3] = 1
                                M_i = np.matrix.getH(ket_i)
                                M_D = np.kron(M_C,M_i)
                                if CHAINLENGTH == 5:
                                    x = x + np.matrix(M_D)*SYSTEM*np.matrix.getH(M_D)
                                if CHAINLENGTH > 5:
                                    for i4 in range(2):
                                        ket_i = np.zeros(2)
                                        ket_i[i4] = 1
                                        M_i = np.matrix.getH(ket_i)
                                        M_E = np.kron(M_D,M_i)
                                        if CHAINLENGTH == 6:
                                            x = x + np.matrix(M_E)*SYSTEM*np.matrix.getH(M_E)
                                        if CHAINLENGTH > 6:
                                            for i5 in range(2):
                                                ket_i = np.zeros(2)
                                                ket_i[i5] = 1
                                                M_i = np.matrix.getH(ket_i)
                                                M_F = np.kron(M_E,M_i)
                                                if CHAINLENGTH == 7:
                                                    x = x + np.matrix(M_F)*SYSTEM*np.matrix.getH(M_F)
                                                if CHAINLENGTH > 7:
                                                    for i6 in range(2):
                                                        ket_i = np.zeros(2)
                                                        ket_i[i6] = 1
                                                        M_i = np.matrix.getH(ket_i)
                                                        M_G = np.kron(M_F,M_i)
                                                        if CHAINLENGTH == 8:
                                                            x = x + np.matrix(M_G)*SYSTEM*np.matrix.getH(M_G)
                                                        if CHAINLENGTH > 8:
                                                            for i7 in range(2):
                                                                ket_i = np.zeros(2)
                                                                ket_i[i7] = 1
                                                                M_i = np.matrix.getH(ket_i)
                                                                M_H = np.kron(M_G,M_i)
                                                                if CHAINLENGTH == 9:
                                                                    x = x + np.matrix(M_H)*SYSTEM*np.matrix.getH(M_H)
                                                                if CHAINLENGTH > 9:
                                                                    for i8 in range(2):
                                                                        ket_i = np.zeros(2)
                                                                        ket_i[i8] = 1
                                                                        M_i = np.matrix.getH(ket_i)
                                                                        M_I = np.kron(M_H,M_i)
                                                                        if CHAINLENGTH == 10:
                                                                            x = x + np.matrix(M_I)*SYSTEM*np.matrix.getH(M_I)
                                                                        if CHAINLENGTH > 10:
                                                                            for i9 in range(2):
                                                                                ket_i = np.zeros(2)
                                                                                ket_i[i9] = 1
                                                                                M_i = np.matrix.getH(ket_i)
                                                                                M_J = np.kron(M_I,M_i)
                                                                                if CHAINLENGTH == 11:
                                                                                    x = x + np.matrix(M_J)*SYSTEM*np.matrix.getH(M_J)           
               
#                x = x + np.matrix(M_)*SYSTEM*np.matrix.getH(M_)
#                y = y + np.matrix(M_C)*SYSTEM1*np.matrix.getH(M_C)
    y = INITIALCHAIN
    z = linalg.sqrtm(np.conj(x-y)*(x-y))
    results.append(0.5*np.absolute(np.trace(z)))
#    z = np.conj(x-y)*(x-y)
    
    SYSTEM = np.matrix(EXP1)*np.matrix(SYSTEM)*np.matrix(EXP2) 
#

#    SYSTEM1 = np.matrix(EXP1)*np.matrix(SYSTEM1)*np.matrix(EXP2)           
    
plt.figure()
plt.plot(results)

