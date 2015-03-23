#NEW QUANTUM WALK/SPIN CHAIN SCRIPT
##################################
#PREAMBLE
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from scipy import linalg
##################################
#INITIAL VARIABLES
STEPS = 100
POSITIONS = 2*STEPS+1
results=[]
results1=[]
I_P = np.identity((POSITIONS))
I_C = np.identity((2))
q0 = np.matrix('1;0')
q1 = np.matrix('0;1')

for x in range(5):
    constant=x/100
    ###################################
    #HADAMARD OPERATOR
    H = (1/np.sqrt(2))*np.matrix([[1,  1],
                                  [1, -1]])
    H = np.kron(I_P,H)
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
    ######################################
    #INITIAL STATE   
    INITIALCOIN =np.matrix([[1],
                            [0]])
    INITIALCOIN = np.outer(INITIALCOIN,INITIALCOIN) 
    
    SYSTEM = np.zeros((2*STEPS+1))
    SYSTEM[STEPS+1] = 1
    SYSTEM = np.outer(SYSTEM,SYSTEM)
    SYSTEM = np.kron(SYSTEM,INITIALCOIN)
    
    ######################################
    #INITIAL STATE
    INITIALCOIN1 =np.matrix([[0],
                             [1]])
    INITIALCOIN1 = np.outer(INITIALCOIN1,INITIALCOIN1) 
    SYSTEM1 = np.zeros((2*STEPS+1))
    SYSTEM1[STEPS+1] = 1
    SYSTEM1 = np.outer(SYSTEM1,SYSTEM1)
    SYSTEM1 = np.kron(SYSTEM1,INITIALCOIN1)
    
    ########################################
    #QUANTUM WALK
    p = np.zeros((POSITIONS,STEPS),dtype=complex)
    results=[]
    for r in range(STEPS):
        
        #Coin Meaurement
        x = np.zeros((2,2))
        y = np.zeros((2,2))
        for i in range(POSITIONS):
            ket_i = np.zeros((POSITIONS,1))
            ket_i[i] = 1
            M_i = np.matrix.getH(ket_i)
            M_P = np.kron(M_i,I_C)
            x = x + np.matrix(M_P)*SYSTEM*np.matrix.getH(M_P)
            y = y + np.matrix(M_P)*SYSTEM1*np.matrix.getH(M_P)
        if r>2:
            z = linalg.sqrtm(np.conj(x-y)*(x-y))
        else:
            z = (np.conj(x-y)*(x-y))        
        results.append(0.5*np.absolute(np.trace(z)))
        
        
    
        E0 = np.matrix([[1, 0], [0, np.sqrt(1-constant)]]) #E1 is the same for both amplitude and phase damping
        E1 = np.matrix([[0, np.sqrt(constant)], [0, 0]])
        E1_PD = np.matrix([[0, 0], [0, np.sqrt(constant)]])
        E0 = np.kron(I_P,E0)
        E1 = np.kron(I_P,E1)
        E1_PD = np.kron(I_P,E1_PD)
    #    system_AD = S*coin_flip*E0*system_AD*np.matrix.getH(E0)*coin_flip*np.matrix.getH(S) + S*coin_flip*E1*system_AD*np.matrix.getH(E1)*coin_flip*np.matrix.getH(S)
    #    system_PD = S*coin_flip*E0*system_PD*np.matrix.getH(E0)*coin_flip*np.matrix.getH(S) + S*coin_flip*E1_PD*system_PD*np.matrix.getH(E1_PD)*coin_flip*np.matrix.getH(S)
    
        SYSTEM = np.matrix(H)*np.matrix(SYSTEM)*np.matrix(np.matrix.getH(H))
        SYSTEM = E0*SYSTEM*np.matrix.getH(E0) + E1*SYSTEM*np.matrix.getH(E1)
        SYSTEM = np.matrix(S)*np.matrix(SYSTEM)*np.matrix(np.matrix.getH(S))  
    #
        SYSTEM1 = np.matrix(H)*np.matrix(SYSTEM1)*np.matrix(np.matrix.getH(H))
        SYSTEM1 = E0*SYSTEM1*np.matrix.getH(E0) + E1*SYSTEM1*np.matrix.getH(E1)
        SYSTEM1 = np.matrix(S)*np.matrix(SYSTEM1)*np.matrix(np.matrix.getH(S))       
    results1.append(results)
     
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

colors = iter(cm.rainbow(np.linspace(0, 1, 5)))
linecycle = ["-","--","-.",":"]
linecycle = iter(linecycle)
y = range(STEPS)
constants = [0,0.01,0.02,0.03,0.04,0.05]

plt.figure(num=1, figsize=(10, 8))
#plt.title('Plot 1', size=14)
plt.xlabel('Position', **font)
plt.ylabel('Probability of Measurement', **font)
for i in range(5):
#    plt.plot(y, x[i], color=next(colors), linestyle=next(linecycle), label=np.around(decoherences_list[i],2))
    plt.plot(y, results1[i], color=next(colors), linestyle="-", label=np.around(constants[i],2))
#plt.plot(y, results, color='#FF0000', linestyle='-', label='Quantum Walk')
#plt.plot(y, x[1], color='#009999', linestyle='-', label='Random Walk')
#plt.plot(x2, y2, color='#FF7400', linestyle=':', label='y2 data')
#plt.plot(x3, y3, color='#00CC00', linestyle='-', label='y2 data')
#plt.plot(x4, y4, color='y', linestyle='-', label='y2 data')
plt.legend(loc='upper right',fancybox=True, fontsize = 'x-small')
#plt.legend(bbox_to_anchor=(0.,1.02,1.,.102), loc = 3,ncol=2,mode='expand', borderaxespad=0.)
#plt.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.,fontsize='small')
plt.show()

