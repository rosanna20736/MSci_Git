from functionheader import quantumwalk
import matplotlib.pyplot as plt
import numpy as np
#Function 'quantumwalk(steps, decoherence type, decoherence value, return type)
#Steps: 50 (Default)
#Decoherence Type: 'n' (none)(Default), 'c' (coin), 'p' (position), 'cp' (coin and position), 'H' (Hadamard), 'm' (Type of measurement chosen at random)
#Decoherence Value: 0 (Default)
#Return Type: 'p' (plot)(Default), 's' (system), 'm' (measurement), 'v' (variance), 'vb' (variance exponent) #'a' (animation)

plt.close('all')
result= []
result1 = []
result2 = []
result3 = []
result4 = []
result5 = []
decoherences_list = []

#INPUT PARAMETERS
decoherences = 18
steps = 100
repeats = 100
repeats1 = 5

#Retrieve variance for coin, position and combined noise types for decoherence prob 0 to 1
for j in range(decoherences+1):
    decoherence = 1 - np.log(1+(np.e-1)*j/decoherences)
    decoherences_list.append(decoherence)
    result.append(quantumwalk(steps,'c', decoherence, 'v'))
    result1.append(quantumwalk(steps,'p', decoherence, 'v'))
    result2.append(quantumwalk(steps,'cp', decoherence, 'v'))
    result3temp = []
    for i in range(repeats):
        result3temp.append(quantumwalk(steps,'H', decoherence, 'v'))
    result3.append(np.average(result3temp))
    result4temp = []
    for i in range(repeats1):
        result4temp.append(quantumwalk(steps,'m', decoherence, 'v'))
    result4.append(np.average(result4temp))
#    result5.append(quantumwalk(steps,'f', decoherence, 'v'))

#plot variance
plt.figure()
plt.plot(decoherences_list,result,label='coin decoherence')
plt.plot(decoherences_list,result1,label='position decoherence')
plt.plot(decoherences_list,result2,label='coin and position decoherence')
plt.plot(decoherences_list,result4,label='random decoherence type',ls='None',marker="o")
plt.plot(decoherences_list,result3,label='imperfect hadamard',ls='None',marker="o")
#plt.plot(decoherences_list,result5,label='analytical approximation')


plt.xlabel('probability of measurement event',fontsize='large')
plt.ylabel('variance of walk',fontsize='large')
#plt.ylabel('exponent',fontsize='large')
plt.legend()
#plt.axis([0,1,0,np.max(result5)])
plt.axis([0,1,1,np.max(result2)])
plt.show

#plot standard deviation
#plt.figure()
#plt.plot(decoherences_list,np.sqrt(result),label='coin decoherence')
#plt.plot(decoherences_list,np.sqrt(result1),label='position decoherence')
#plt.plot(decoherences_list,np.sqrt(result2),label='coin and position decoherence')
#plt.plot(decoherences_list,np.sqrt(result4),label='random decoherence type')
#plt.plot(decoherences_list,np.sqrt(result3),label='imperfect hadamard',ls='None',marker="o")

#plt.xlabel('probability of measurement event',fontsize='large')
#plt.ylabel('standard deviation of walk',fontsize='large')
#plt.legend()
#plt.show
#plt.axis([0,1,0,np.max(np.sqrt(result))])

#for plotting several distributions on 1 plot - replace 'v' with 'p' in one of the  functions above
#plt.xlabel('particle position, x',fontsize='large')
#plt.ylabel('probability particle is found at position x',fontsize='large')
#plt.legend()
#plt.show