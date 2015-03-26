#NEW QUANTUM WALK/SPIN CHAIN SCRIPT
##################################
#PREAMBLE
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib
from scipy import linalg
##################################
#INITIAL VARIABLES
STEPS = 100
CHAINLENGTH = 3
results=[]
J = 0.1
J1 = 0.05
I_C = np.identity((2))
chaintype = 1
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
EXP1 = linalg.expm2(-1*A)
EXP2 = linalg.expm2(A)
######################################
INITIALCOIN = np.matrix([[1],
                         [0]])
INITIALCOIN = np.outer(INITIALCOIN,INITIALCOIN)
INITIALCHAIN = np.matrix([[0],
                          [1]])
#
INITIALCHAIN = np.outer(INITIALCHAIN,INITIALCHAIN)
A = np.kron(INITIALCOIN,INITIALCHAIN)                                         
for X in range(CHAINLENGTH-2):
    A = np.kron(A,INITIALCHAIN)
SYSTEM = A
########################################
INITIALCOIN2 = np.matrix([[0],
                         [1]])
INITIALCOIN2 = np.outer(INITIALCOIN2,INITIALCOIN2)
INITIALCHAIN2 = np.matrix([[0],
                          [1]])
INITIALCHAIN2 = np.outer(INITIALCHAIN2,INITIALCHAIN2)
A2 = np.kron(INITIALCOIN2,INITIALCHAIN2)                                         
for Y in range(CHAINLENGTH-2):
    A2 = np.kron(A2,INITIALCHAIN2)
SYSTEM2 = A2
########################################
#QUANTUM WALK
ket1 = np.matrix([[1],[0]])
ket2 = np.matrix([[0],[1]])
for r in range(STEPS):
    NEW1 = SYSTEM
    NEW2 = SYSTEM2
    for X in range(CHAINLENGTH,1,-1):
        OPP = I_C
        NEWSYSTEM = NEW1
        for Y in range(X-2):
            OPP = np.kron(OPP,I_C)
        OPP1 = np.kron(OPP,ket1)
        OPP2 = np.kron(OPP,ket2)
        OPP3 = np.kron(OPP,np.matrix.getH(ket1))
        OPP4 = np.kron(OPP,np.matrix.getH(ket2))
        dim = (np.power(2,X-1),np.power(2,X-1))
        NEW1 = np.zeros(dim)
        NEW1 = OPP3*NEWSYSTEM*OPP1 + OPP4*NEWSYSTEM*OPP2
    for v in range(CHAINLENGTH,1,-1):     
        OPP = I_C
        NEWSYSTEM2 = NEW2
        for w in range(v-2):
            OPP = np.kron(OPP,I_C)
        OPP1 = np.kron(OPP,ket1)
        OPP2 = np.kron(OPP,ket2)
        OPP3 = np.kron(OPP,np.matrix.getH(ket1))
        OPP4 = np.kron(OPP,np.matrix.getH(ket2))
        dim = (np.power(2,v-1),np.power(2,v-1))
        NEW2 = np.zeros(dim)
        NEW2 = OPP3*NEWSYSTEM2*OPP1 + OPP4*NEWSYSTEM2*OPP2
    x= NEW1  
    y= NEW2
#    y = INITIALCHAIN
    z = linalg.sqrtm(np.matrix.getH(x-y)*(x-y))
    results.append(0.5*np.absolute(np.trace(z)))
    
    SYSTEM = np.matrix(EXP1)*np.matrix(SYSTEM)*np.matrix(EXP2)            
    SYSTEM2 = np.matrix(EXP1)*np.matrix(SYSTEM2)*np.matrix(EXP2)   
###################################PLOT OPTIONS    

font = {'family' : 'sans-serif',
'weight' : 'light',
'size' : 12}
ax = plt.subplot(111)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.tick_params(axis="both", which="both", bottom="on", top="off",
labelbottom="on", left="on", right="off", labelleft="on")
ax.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
ax.grid(axis='x', color="0.9", linestyle='-', linewidth=1)
#ax.yaxis.grid(True)
#ax.xaxis.grid(True)

matplotlib.rc('font', **font)

ax.set_axisbelow(True) 
plt.figure(num=1, figsize=(10, 8))
plt.ylim(0,1)
plt.xlim(0,STEPS-1)
#plt.title('Plot 1', size=14)
plt.xlabel('Steps', **font)
plt.ylabel('Trace Distance', **font)
plt.plot(results, color='#FF0000', linestyle='-', label= 'TD J(%s) J1(%s)'%(J, J1))
#plt.plot(x1, y1, color='#FF7400', linestyle='--', label='y2 data')
#plt.plot(x2, y2, color='#009999', linestyle=':', label='y2 data')
#plt.plot(x3, y3, color='#00CC00', linestyle='-', label='y2 data')
#plt.plot(x4, y4, color='y', linestyle='-', label='y2 data')
plt.legend(loc='lower right',fancybox=True)

plt.show()


