#NEW QUANTUM WALK/SPIN CHAIN SCRIPT
##################################
#PREAMBLE
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from spinchaingenerator import U_spinchain
import time
##################################
print('running')
plt.close('all')
#INITIAL VARIABLES
STEPS = 500
CHAINLENGTH = 5
results=[]
J=0.2
J1=0.1
I_C = np.identity((2))
chaintype = 0
##################################
#SPINCHAIN OPERATOR
pauli_x = np.matrix(('0 1;1 0'))
pauli_y = np.matrix(('0 -1j;1j 0'))
pauli_z = np.matrix(('1 0;0 -1'))
pauli_i = np.matrix(('1 0;0 1'))
if chaintype == 1:
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
    EXP1 = linalg.expm2(-1*A)
    EXP2 = linalg.expm2(A)
if chaintype == 0:
    EXP2 = U_spinchain('xy',CHAINLENGTH,J,J1)
    EXP1 = np.matrix.getH(EXP2)
if chaintype == 2:
    A = (np.kron(np.kron(np.kron(np.kron(pauli_x,pauli_x),I_C),I_C),I_C)*J1)+ (np.kron(np.kron(np.kron(np.kron(I_C,pauli_x),pauli_x),I_C),I_C)*J)+ (np.kron(np.kron(np.kron(np.kron(I_C,I_C),pauli_x),pauli_x),I_C)*J)+ (np.kron(np.kron(np.kron(np.kron(I_C,I_C),I_C),pauli_x),pauli_x)*J)+ (np.kron(np.kron(np.kron(np.kron(pauli_y,pauli_y),I_C),I_C),I_C)*J1)+ (np.kron(np.kron(np.kron(np.kron(I_C,pauli_y),pauli_y),I_C),I_C)*J)+ (np.kron(np.kron(np.kron(np.kron(I_C,I_C),pauli_y),pauli_y),I_C)*J)+ (np.kron(np.kron(np.kron(np.kron(I_C,I_C),I_C),pauli_y),pauli_y)*J)
    A = np.matrix(A)*1j
    EXP1 = linalg.expm2(-1*A)
    EXP2 = linalg.expm2(A)

######################################
INITIALCOIN = np.matrix([[0],
                         [1]])
INITIALCOIN_orth = np.matrix([[1],
                         [0]])
INITIALCOIN = np.outer(INITIALCOIN,INITIALCOIN)
INITIALCOIN_orth = np.outer(INITIALCOIN_orth,INITIALCOIN_orth)
INITIALCHAIN = np.matrix([[1],
                          [0]])
INITIALCHAIN = np.outer(INITIALCHAIN,INITIALCHAIN)
A = np.kron(INITIALCOIN,INITIALCHAIN)
A_orth = np.kron(INITIALCOIN_orth,INITIALCHAIN)                                       
for x in range(CHAINLENGTH-2):
    A = np.kron(A,INITIALCHAIN)
    A_orth = np.kron(A_orth,INITIALCHAIN)
SYSTEM = A
SYSTEM_orth = A_orth
########################################
start_time = time.time()
#QUANTUM WALK
ket1 = np.matrix([[1],[0]])
ket2 = np.matrix([[0],[1]])
for r in range(STEPS):
    NEW1 = SYSTEM
    for x in range(CHAINLENGTH,1,-1):
        OPP = I_C
        NEWSYSTEM = NEW1
        for y in range(x-2):
            OPP = np.kron(OPP,I_C)
        OPP1 = np.kron(OPP,ket1)
        OPP2 = np.kron(OPP,ket2)
        OPP3 = np.kron(OPP,np.matrix.getH(ket1))
        OPP4 = np.kron(OPP,np.matrix.getH(ket2))
        dim = (np.power(2,x-1),np.power(2,x-1))
        NEW1 = np.zeros(dim)
        NEW1 = OPP3*NEWSYSTEM*OPP1 + OPP4*NEWSYSTEM*OPP2
    x = NEW1
    y = INITIALCHAIN
    z = linalg.sqrtm(np.conj(x-y)*(x-y))
    results.append(0.5*np.absolute(np.trace(z)))
    
    SYSTEM = np.matrix(EXP1)*np.matrix(SYSTEM)*np.matrix(EXP2)          
print("Michael's =  %s seconds" % (time.time() - start_time))
plt.figure()
plt.plot(results)

SYSTEM = A
results = []
start_time = time.time()
ket0 = np.matrix([[1],[0]])
ket1 = np.matrix([[0],[1]])
rho_orth = np.kron(ket1,np.matrix.getH(ket1))
for r in range(STEPS):
    rho_C = np.zeros((2,2))
    rho_orth = np.zeros((2,2))
    for i in range(np.power(2,CHAINLENGTH-1)):
        ket_i = np.zeros((np.power(2,CHAINLENGTH-1),1))
        ket_i[i] = 1
        OPP_left = np.kron(I_C,np.matrix.getH(ket_i))
        OPP_right = np.kron(I_C,ket_i)
        rho_C = rho_C + np.matrix(OPP_left)*np.matrix(SYSTEM)*np.matrix(OPP_right)
        rho_orth = rho_orth + np.matrix(OPP_left)*np.matrix(SYSTEM_orth)*np.matrix(OPP_right)
    rho_diff = np.matrix(rho_C)-np.matrix(rho_orth)
    results.append(0.5*np.trace(linalg.sqrtm(np.matrix.getH(rho_diff)*np.matrix(rho_diff))))
    SYSTEM = np.matrix(EXP1)*np.matrix(SYSTEM)*np.matrix(EXP2)   
    SYSTEM_orth = np.matrix(EXP1)*np.matrix(SYSTEM_orth)*np.matrix(EXP2)  
print("Rosanna's =  %s seconds" % (time.time() - start_time))
plt.figure()
plt.plot(results)

