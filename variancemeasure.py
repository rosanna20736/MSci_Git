from functionheader import quantumwalk
import matplotlib.pyplot as plt
import numpy as np
#Function 'quantumwalk(steps, decoherence type, decoherence value, return type)
#Steps: 50 (Default)
#Decoherence Type: 'n' (none)(Default), 'c' (coin), 'p' (position), 'cp' (coin and position), 'H' (Hadamard) 
#Decoherence Value: 0 (Default)
#Return Type: 'p' (plot)(Default), 's' (system), 'm' (measurement), 'v' (variance), #'a' (animation)

plt.close('all')
result= []
result1 = []
result2 = []
decoherences_list = []

#INPUT PARAMETERS
decoherences = 50
steps = 50

#Retrieve variance for coin, position and combined noise types for decoherence prob 0 to 1
for j in range(decoherences+1):
    decoherence = 1 - np.log(1+(np.e-1)*j/decoherences)
    decoherences_list.append(decoherence)
    result.append(quantumwalk(steps,'c', decoherence, 'v'))
    result1.append(quantumwalk(steps,'p', decoherence, 'v'))
    result2.append(quantumwalk(steps,'cp', decoherence, 'v'))

#plot variance
plt.figure(1)
plt.plot(decoherences_list,result,label='coin decoherence')
plt.plot(decoherences_list,result1,label='position decoherence')
plt.plot(decoherences_list,result2,label='coin and position decoherence')

plt.xlabel('probability of measurement event',fontsize='large')
plt.ylabel('variance of walk',fontsize='large')
plt.legend()
plt.show

#plot standard deviation
plt.figure(2)
plt.plot(decoherences_list,np.sqrt(result),label='coin decoherence')
plt.plot(decoherences_list,np.sqrt(result1),label='position decoherence')
plt.plot(decoherences_list,np.sqrt(result2),label='coin and position decoherence')

plt.xlabel('probability of measurement event',fontsize='large')
plt.ylabel('standard deviation of walk',fontsize='large')
plt.legend()
plt.show