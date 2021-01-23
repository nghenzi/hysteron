# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 13:35:18 2016

@author: marian
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

plt.close('all')

filename="icc_100nA.dat"
data = np.loadtxt(filename, delimiter="\t", skiprows=1, usecols=[0,1])

file_1nA="non_volatile.dat"
data_1nA = np.loadtxt(file_1nA, delimiter="\t", skiprows=1, usecols=[0,1])


current_data = []
voltage_data = []
current_data_1nA = []
voltage_data_1nA = []

voltage_data= data[:,0]     
current_data= data[:,1]
voltage_data_1nA= data_1nA[:,0]     
current_data_1nA= data_1nA[:,1]

# Voltage driven memristor
# current at READ with 0.1V (in A)
i0_max_pos = 4e-7
i0_max_neg = 4e-7
i0_max = 1
i0_min = 1e-15

# Voltage thresholds for switching
v_threshold_pos = 5.5
v_off_pos=.1
v_threshold_neg = -7.5
v_off_neg = -.1

# some constants
r_pos = 5.0
r_neg = 3.0
R1 = float('inf')
R2 = float('inf')
G1 = 1.0/R1
G2 = 1.0/R2
alpha_pos =2.1
alpha_neg=1.1
alpha=1  
Rs = 29.0e6
d = 1.0/(1 + Rs*G1)
beta_th=1.0



#T = 2*np.pi
#t = np.linspace(0, 2*T, 100)
#t0 = 0.5*np.pi
#v_sweep = 2*sig.sawtooth(t + t0, 0.5)

# OJO QUE HAY UN PROBLEMA DE CONVERGENCIA CON AMPLITUDES ALTAS EN EL NEWTON RAPHSON
# EL PROBLEMA ES UNA RELACIÓN ENTRE EL ALPHA, LA AMPLITUD DE LA SEÑAL EN LA EXPONENCIAL



# voltage input signal
MAX = 20.0
MIN = -20.0

pos_up = np.linspace(0, MAX, 200)
pos_down = pos_up[::-1]
pos = np.concatenate((pos_up, pos_down))
neg_down = np.linspace(0, MIN, 200)
neg_up = neg_down[::-1]
neg = np.concatenate((neg_down, neg_up))

v_sweep = np.concatenate((pos, neg))
#graph = plt.figure(0)
#plt.plot(v_sweep, 'bo')
#plt.grid(True)




## 100 nA   lamb_max=1.0   y beta_th=1
## 10 nA    lamb_max= .27   y beta_th=.74
## 1nA      lamb_max= .07   y beta_th=.55
## 100 pA   lamb_max= .008   y beta_th= .5

def hysteron(v_sweep,Lmax):
    # allocating
    lamb_array = []
    i0_array = []
    w_array = []
    current_array = []
    voltage_array = []
    beta_th=1*Lmax    
    # initial condicionts for flags (to indicate if MAX or MIN are reached)
    flag_max = False
    flag_min = False
    for v in v_sweep:
        if v >= 0:
            alpha=alpha_pos
            i0_max=i0_max_pos*Lmax
            if flag_max and not flag_min:
                lamb = Lmax*1.0/(1.0 + np.exp(-r_pos*(v - v_off_pos/beta_th)))
            else:
                lamb = Lmax*1.0/(1.0 + np.exp(-r_pos*(v - v_threshold_pos))**3)
            if v == MAX:
                flag_max = True
                flag_min = False
        elif v < 0:
            alpha=alpha_neg
            i0_max=i0_max_neg*Lmax
            if flag_max and not flag_min:
                lamb = Lmax*1.0/(1.0 + np.exp(r_neg*(v - v_threshold_neg)))
            else:
                lamb = Lmax*1.0/(1.0 + np.exp(r_neg*(v - v_off_neg/beta_th)))
            if v == MIN:
                flag_max = False
                flag_min = True
        i0 = i0_min + (i0_max - i0_min)*lamb*.01
        # Rs and alpha are constants in this approach
        jk = alpha*i0*Rs*d*np.exp(alpha*d*(abs(v) + Rs*i0))
        w = np.log(1+jk)*( 1- ( (np.log(1+(np.log(1+jk))))/(2+(np.log(1+jk))) ) )
        
        current =np.absolute( np.sign(v)*(w/(alpha*Rs) + d*(G1*abs(v)-i0) + G2*abs(v)))+abs(v)/2e13+6e-13
#        current1 =np.absolute( np.sign(v)*(w/(alpha*Rs) + d*(G1*abs(v)-i0) + G2*abs(v))+abs(v)/2e13+6e-13)
        # storing in arrays    
        lamb_array.append(lamb)
        i0_array.append(i0)
        w_array.append(w)
        current_array.append(current)
#        if current1> 5e1:
#            current1=5e-11
#        current_array1.append(current1)
        voltage_array.append(v)
    return voltage_array, current_array,lamb_array


v1, i1, l1 = hysteron(v_sweep,1)

v2, i2, l2 = hysteron(v_sweep,0.1)  

v3, i3, l3 = hysteron(v_sweep,0.05)  
# plot all data
graph = plt.figure(11)
#plt.semilogy(voltage_array, current_array,  'bo-', label = 'original')
#plt.plot(voltage_array, current_array1, '-', label = 'simulation')

#plt.semilogy(voltage_data, current_data, 'go', label = 'exp_100nA')
plt.plot(v1, l1, 'r-', linewidth=3.0,  label = 'simulation')
#plt.plot(voltage_data_1nA, current_data_1nA, 'bo-.', label = 'exp_1nA')
plt.plot(v2, l2, 'g-', linewidth=3.0,  label = 'simulation')
plt.plot(v3, l3, 'y-', linewidth=3.0,  label = 'simulation')
plt.legend(loc = 'best')
plt.grid(True)
plt.show()

graph = plt.figure(1)
#plt.semilogy(voltage_array, current_array,  'bo-', label = 'original')
plt.semilogy(v1, i1, 'r-', linewidth=3.0,  label = 'simulation')
plt.semilogy(v2, i2, 'g-', linewidth=3.0,  label = 'simulation')
plt.plot(v3, i3, 'y-', linewidth=3.0,  label = 'simulation')

plt.semilogy(voltage_data_1nA, current_data_1nA, 'bo-.', label = 'exp_1nA')
plt.legend(loc = 'best')
plt.grid(True)
plt.show()



f = open ("fichero0553.dat", "w")


for i in range(len(v1)):
    f.write("%5.18f %5.18f\n" % (v1[i], i1[i]))
f.close()

f = open ("fichero0554.dat", "w")


for i in range(len(v1)):
    f.write("%5.18f %5.18f\n" % (v2[i], i2[i]))
f.close()

f = open ("fichero0555.dat", "w")


for i in range(len(v1)):
    f.write("%5.18f %5.18f\n" % (v3[i], i3[i]))
f.close()

print ('End of simulation')