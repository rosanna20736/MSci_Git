#NEW QUANTUM WALK/SPIN CHAIN SCRIPT
##################################
#PREAMBLE
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import linalg
import time
##################################
#INITIAL VARIABLES
STEPS = 100
POSITIONS = 2*STEPS+1
CHAINLENGTH = 2
results=[]
times = []
J = 0.1
J1 =0.1
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
for z in range(CHAINLENGTH-1):
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
A = (J1*CHAIN1 + J*CHAIN2)
A = 1j*A
EXP1 = linalg.expm2(-1*A)
EXP2 = linalg.expm2(A)
EXP1 = np.kron(I_P,EXP1)
EXP2 = np.kron(I_P,EXP2)
######################################
#INITIAL STATE   
INITIALCOIN =               (1/np.sqrt(2))*np.matrix([[1],
                                       [1]])
INITIALCOIN = np.outer(INITIALCOIN,INITIALCOIN)
INITIALCHAIN = np.matrix([[1,  0],
                          [0,  0]])  
SYSTEM = np.zeros((2*STEPS+1))
SYSTEM[STEPS+1] = 1
SYSTEM = np.outer(SYSTEM,SYSTEM)
A = np.kron(INITIALCOIN,INITIALCHAIN)
for z in range(CHAINLENGTH-2):
    A = np.kron(A,INITIALCHAIN)
SYSTEM = np.kron(SYSTEM,A)

######################################
#INITIAL STATE
INITIALCOIN1 =               (1/np.sqrt(2))*np.matrix([[1],
                                       [-1]])
INITIALCOIN1 = np.outer(INITIALCOIN1,INITIALCOIN1)
INITIALCHAIN1 = np.matrix([[1,  0],
                           [0,  0]])  
SYSTEM1 = np.zeros((2*STEPS+1))
SYSTEM1[STEPS+1] = 1
SYSTEM1 = np.outer(SYSTEM1,SYSTEM1)
A1 = np.kron(INITIALCOIN1,INITIALCHAIN1)
for z in range(CHAINLENGTH-2):
    A1 = np.kron(A1,INITIALCHAIN1)
SYSTEM1 = np.kron(SYSTEM1,A1)
########################################
#QUANTUM WALK

ket1 = np.matrix([[1],[0]])
ket2 = np.matrix([[0],[1]])
for r in range(STEPS):
    start = time.time()
    NEW1 = SYSTEM
    NEW2 = SYSTEM1
    for X in range(CHAINLENGTH,1,-1):
        OPP = np.kron(I_P,I_C)
        NEWSYSTEM = NEW1
        for Y in range(X-2):
            OPP = np.kron(OPP,I_C)
        OPP1 = np.kron(OPP,ket1)
        OPP2 = np.kron(OPP,ket2)
        OPP3 = np.kron(OPP,np.matrix.getH(ket1))
        OPP4 = np.kron(OPP,np.matrix.getH(ket2))
        dim = (np.power(2,X-1)*POSITIONS,np.power(2,X-1)*POSITIONS)
        NEW1 = np.zeros(dim)
        NEW1 = OPP3*NEWSYSTEM*OPP1 + OPP4*NEWSYSTEM*OPP2
    for v in range(CHAINLENGTH,1,-1):     
        OPP = np.kron(I_P,I_C)
        NEWSYSTEM2 = NEW2
        for w in range(v-2):
            OPP = np.kron(OPP,I_C)
        OPP1 = np.kron(OPP,ket1)
        OPP2 = np.kron(OPP,ket2)
        OPP3 = np.kron(OPP,np.matrix.getH(ket1))
        OPP4 = np.kron(OPP,np.matrix.getH(ket2))
        dim = (np.power(2,v-1)*POSITIONS,np.power(2,v-1)*POSITIONS)
        NEW2 = np.zeros(dim)
        NEW2 = OPP3*NEWSYSTEM2*OPP1 + OPP4*NEWSYSTEM2*OPP2
    x= NEW1  
    y= NEW2
#    y = INITIALCHAIN
#    z = linalg.sqrtm(np.matrix.getH(x-y)*(x-y))
    z = (np.matrix.getH(x-y)*(x-y))
    results.append(0.5*np.absolute(np.trace(z)))
    SYSTEM = np.matrix(H)*np.matrix(SYSTEM)*np.matrix(np.matrix.getH(H))
    SYSTEM = np.matrix(S)*np.matrix(SYSTEM)*np.matrix(np.matrix.getH(S))
    SYSTEM = np.matrix(EXP2)*np.matrix(SYSTEM)*np.matrix(EXP1)     
#
    SYSTEM1 = np.matrix(H)*np.matrix(SYSTEM1)*np.matrix(np.matrix.getH(H))
    SYSTEM1 = np.matrix(S)*np.matrix(SYSTEM1)*np.matrix(np.matrix.getH(S))
    SYSTEM1 = np.matrix(EXP1)*np.matrix(SYSTEM1)*np.matrix(EXP2) 
          
    end = time.time()
    times.append(end-start)
    average = sum(times)/len(times)
    timeleft = (STEPS - r)*average
    print((sum(times)/(sum(times)+timeleft))*100)
    
     
######################################PLOTOPTIONS###############################

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
ax.set_axisbelow(True)
matplotlib.rc('font', **font)

#colors = iter(cm.rainbow(np.linspace(0, 1, decoherences+1)))
linecycle = ["-","--","-.",":"]
linecycle = iter(linecycle)
y = range(STEPS)

plt.figure(num=1, figsize=(10, 8))
#plt.title('Plot 1', size=14)
plt.xlabel('steps', **font)
plt.ylabel('Trace Distance', **font)
#for i in range(decoherences+1):
#    plt.plot(y, x[i], color=next(colors), linestyle=next(linecycle), label=np.around(decoherences_list[i],2))
#    plt.plot(y, x[i], color=next(colors), linestyle="-", label=np.around(decoherences_list[i],2))
plt.plot(y, results, color='#FF0000', linestyle='-', label='Quantum Walk')
#plt.plot(y, x[1], color='#009999', linestyle='-', label='Random Walk')
#plt.plot(x2, y2, color='#FF7400', linestyle=':', label='y2 data')
#plt.plot(x3, y3, color='#00CC00', linestyle='-', label='y2 data')
#plt.plot(x4, y4, color='y', linestyle='-', label='y2 data')
plt.legend(loc='upper right',fancybox=True, fontsize = 'x-small')
#plt.legend(bbox_to_anchor=(0.,1.02,1.,.102), loc = 3,ncol=2,mode='expand', borderaxespad=0.)
#plt.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.,fontsize='small')
plt.show()


