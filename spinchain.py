#import libraries
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.animation as animation
from scipy import linalg

chainlength = 3
JT = 1

pauli_x = np.matrix('0 1;1 0')
pauli_y = np.matrix('0 -1j;1j 0')
pauli_z = np.matrix('1 0;0 -1')
pauli_i = np.matrix('1 0;0 1')

def sigma_n(pauli,n):   #LABELS ON SPIN CHAIN GO FROM particle 0 TO particle N-1
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

def U_spinchain(xyz):
    #function called to get exp(-iH) where H is the spin chain hamiltonian
    H = 0
    expH = 0
    if xyz == 'z':
        for i in range(chainlength-1):
            H = H + sigma_n('z',i)*sigma_n('z',i+1)
        H = H + sigma_n('z',chainlength-1)*sigma_n('z',0)
    elif xyz == 'xy':
        for i in range(chainlength-1):
            H = H + sigma_n('x',i)*sigma_n('x',i+1) + sigma_n('y',i)*sigma_n('y',i+1)
        H = H + sigma_n('x',chainlength-1)*sigma_n('x',0) + sigma_n('y',chainlength-1)*sigma_n('y',0)
    else:
        print('H_spinchain: invalid input')
    H = H*JT
    expH = linalg.expm2(-1j*H)
    return expH
    
#test they are unitary i.e. test1 and test2 should be identity
test1 = U_spinchain('z')*np.matrix.getH(U_spinchain('z'))
test2 = U_spinchain('xy')*np.matrix.getH(U_spinchain('xy'))