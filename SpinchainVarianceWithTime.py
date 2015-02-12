#NEW QUANTUM WALK/SPIN CHAIN SCRIPT
##################################
#PREAMBLE
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
from scipy import linalg
##################################
print('running')
plt.close('all')
#INITIAL VARIABLES
STEPS = 50
POSITIONS = 2*STEPS+1
CHAINLENGTH = 2
results=[]
results_Q=[]
J = 0.1
J1 = 0.1
I_P = np.identity((POSITIONS))
I_C = np.identity((2))
##################################
#ANIMATION SET-UP
#fig = plt.figure()
#ax = plt.axes(xlim=(-STEPS, STEPS), ylim=(0, 0.5))
#line, = ax.plot([], [], lw=2)
#
#
#def init():
#    line.set_data([], [])
#    return line,
##################################
#HADAMARD OPERATOR
H = (1/np.sqrt(2))*np.matrix([[1,  1],
                              [1, -1]])
H = np.kron(I_P,H)
H_nochain = H
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

SHIFT_nochain = SHIFTUP+SHIFTDOWN

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
SYSTEM = np.kron(SYSTEM,I_C)
SYSTEM_nochain = SYSTEM
for ab in range(CHAINLENGTH-1):
    SYSTEM = np.kron(SYSTEM,I_C)
    
INITIALCOIN =               np.matrix([[1],
                                       [0]])
INITIALCOIN = np.outer(INITIALCOIN,INITIALCOIN)
INITIALCOIN = np.kron(I_P,INITIALCOIN)
INITIALCOIN_nochain = INITIALCOIN
for aC in range(CHAINLENGTH-1):
    INITIALCOIN = np.kron(INITIALCOIN,I_C)    

INITIALCHAIN = np.matrix([[0],
                          [1]])
INITIALCHAIN1 = np.outer(INITIALCHAIN,INITIALCHAIN)
INITIALCHAIN = np.kron(I_P,I_C)        
for ad in range(CHAINLENGTH-1):
    INITIALCHAIN = np.kron(INITIALCHAIN,INITIALCHAIN1)   

SYSTEM = np.matrix(SYSTEM)*np.matrix(INITIALCOIN)*np.matrix(INITIALCHAIN)
SYSTEM_Q = np.matrix(SYSTEM_nochain)*np.matrix(INITIALCOIN_nochain)
#print(SYSTEM_C.shape)
del INITIALCHAIN
del INITIALCOIN
del INITIALCOIN_nochain
del INITIALCHAIN1
##################################
############################################################ 
#LAST BIT OF CHAIN CODE
CHAIN = J1*CHAIN1 + J*CHAIN2
CHAIN = np.kron(I_P,CHAIN)
EXP1 = linalg.expm2(-1j*CHAIN)
EXP2 = linalg.expm2(1j*CHAIN)
#EXP1 = np.kron(I_P,EXP1)
#EXP2 = np.kron(I_P,EXP2)
###############################################################
#QUANTUM WALK
p = np.zeros((POSITIONS,STEPS),dtype=complex)
p_Q = np.zeros((POSITIONS,STEPS),dtype=complex)
for r in range(STEPS):
    SYSTEM = np.matrix(EXP1)*np.matrix(SYSTEM)*np.matrix(EXP2)    
    SYSTEM = np.matrix(H)*np.matrix(SYSTEM)*np.matrix(np.matrix.getH(H))
    SYSTEM = np.matrix(SHIFT)*np.matrix(SYSTEM)*np.matrix.getH(SHIFT)
    SYSTEM_Q = np.matrix(H_nochain)*np.matrix(SYSTEM_Q)*np.matrix.getH(H_nochain)
    SYSTEM_Q = np.matrix(SHIFT_nochain)*np.matrix(SYSTEM_Q)*np.matrix.getH(SHIFT_nochain)
##################################
    #MEASUREMENT
    
    for i in range(POSITIONS):
        ket_i = np.zeros((POSITIONS,1))
        ket_i[i] = 1
        M_i = np.outer(ket_i,ket_i)
        M_P = np.kron(M_i,I_C)
        M_P_nochain = M_P
        for q in range(CHAINLENGTH-1):
            M_P = np.kron(M_P,I_C)
        p[i,r]= np.trace(M_P*SYSTEM*np.matrix.getH(M_P))
        p_Q[i,r]= np.trace(M_P_nochain*SYSTEM_Q*np.matrix.getH(M_P_nochain))
#    if r%2==0:################################################################CODETOBEFIXED###########################
#        COPY = list(np.absolute(p[0:POSITIONS:2,r]))
#    else:
#        COPY = list(np.absolute(p[1:POSITIONS:2,r]))
#    LABELS=list(range(-STEPS,STEPS+2,2))
#    ##################################
    COPY = list(np.absolute(p[:,r]))
    COPY_Q = list(np.absolute(p_Q[:,r]))
    LABELS = list(range(-STEPS,STEPS))
    #ANALYSIS
    variance = 0
    mean = 0
    variance_Q = 0
    mean_Q = 0
    for i in range(0,2*STEPS):
        mean = mean + COPY[i]*LABELS[i]
        mean_Q = mean_Q + COPY_Q[i]*LABELS[i]
    mean = mean/(2*r+1)
    mean_Q = mean_Q/(2*r+1)
    for j in range(0,2*STEPS):
        variance = variance + COPY[j]*(LABELS[j] - mean)*(LABELS[j] - mean)
        variance_Q = variance_Q + COPY_Q[j]*(LABELS[j] - mean_Q)*(LABELS[j] - mean_Q)
    results.append(variance)
    results_Q.append(variance_Q)
##################################
    
#def animate(r): 
#    line.set_data(LABELS, p[0:2*STEPS+1:2,r])
#    return line,
#    
#anim = animation.FuncAnimation(fig, animate, init_func=init, frames=STEPS, interval=50, blit=True)
#plt.show()
    
plt.plot(results_Q,label='fully quantum walk')
plt.plot(results,label='J = ' + str(round(J,2)) + ', chain length = ' + str(CHAINLENGTH))
plt.legend(loc=2)
plt.xlabel('number of steps',fontsize='large')
plt.ylabel('variance of walk',fontsize='large')
plt.show