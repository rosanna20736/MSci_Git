#import libraries
import numpy as np
from scipy import linalg
print('running spinchaingenerator.py')

pauli_x = np.matrix('0 1;1 0')
pauli_y = np.matrix('0 -1j;1j 0')
pauli_z = np.matrix('1 0;0 -1')
pauli_i = np.matrix('1 0;0 1')

def sigma_n(pauli,n,chainlength):   #LABELS ON SPIN CHAIN GO FROM particle 0 TO particle N-1
    #function acting pauli x, y or z on nth particle and identity on others
    sigma = 1
    if pauli == 'x':
        for i in range(chainlength):
            if i==n:
                sigma = np.kron(sigma,pauli_x)
            else:
                sigma = np.kron(sigma,pauli_i)
    elif pauli == 'y':
        for i in range(chainlength):
            if i==n:
                sigma = np.kron(sigma,pauli_y)
            else:
                sigma = np.kron(sigma,pauli_i)
    elif pauli == 'z':
        for i in range(chainlength):
            if i==n:
                sigma = np.kron(sigma,pauli_z)
            else:
                sigma = np.kron(sigma,pauli_i)
    else:
        print('sigma_n: invalid input')
    return sigma

def U_spinchain(xyz,chainlength,J,J1):
    #function called to get exp(-iH) where H is the spin chain hamiltonian
    H = 0
    expH = 0
    if xyz == 'z':
        H = J1*sigma_n('z',0,chainlength)*sigma_n('z',1,chainlength)
        for i in range(1,chainlength-1):
            H = H + J*sigma_n('z',i,chainlength)*sigma_n('z',i+1,chainlength)
    elif xyz == 'xy':
        H = J1*(sigma_n('x',0,chainlength)*sigma_n('x',1,chainlength) + sigma_n('y',0,chainlength)*sigma_n('y',1,chainlength))
        for i in range(1,chainlength-1):
            H = H + J*(sigma_n('x',i,chainlength)*sigma_n('x',i+1,chainlength) + sigma_n('y',i,chainlength)*sigma_n('y',i+1,chainlength))
        #for closed system where nth qubit interacts with first qubit
        #H = H + J*(sigma_n('x',chainlength-1,chainlength)*sigma_n('x',0,chainlength) + sigma_n('y',chainlength-1,chainlength)*sigma_n('y',0,chainlength))
    else:
        print('H_spinchain: invalid input')
    expH = linalg.expm2(-1j*H)
    return expH