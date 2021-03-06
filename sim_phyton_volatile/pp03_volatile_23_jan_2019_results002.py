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

# Load data from tsv file: data
data = np.loadtxt(filename, delimiter="\t", skiprows=1, 
                  usecols=[0,1])


current_data = []
voltage_data = []

voltage_data= data[:,0]     
current_data= data[:,1]

# Voltage driven memristor
# current at READ with 0.1V (in A)
i0_max = 0.00000001
i0_min = 0.000000000000001

# Voltage thresholds for switching
v_threshold_pos = 4.6
v_off_pos=1.3
v_threshold_neg = -2.0
v_off_neg = -1.0

# some constants
r_pos = 25.0
r_neg = 25.0
R1 = float('inf')
R2 = 1e14
G1 = 1.0/R1
G2 = 1.0/R2
alpha =2  
Rs = 45000000.0
d = 1.0/(1 + Rs*G1)

# allocating
lamb_array = []
i0_array = []
w_array = []
current_array = []
current_array1 = []
voltage_array = []

#T = 2*np.pi
#t = np.linspace(0, 2*T, 100)
#t0 = 0.5*np.pi
#v_sweep = 2*sig.sawtooth(t + t0, 0.5)

# OJO QUE HAY UN PROBLEMA DE CONVERGENCIA CON AMPLITUDES ALTAS EN EL NEWTON RAPHSON
# EL PROBLEMA ES UNA RELACIÓN ENTRE EL ALPHA, LA AMPLITUD DE LA SEÑAL EN LA EXPONENCIAL

# voltage input signal
MAX = 10.0
MIN = 0
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

# initial condicionts for flags (to indicate if MAX or MIN are reached)
flag_max = False
flag_min = False

for v in v_sweep:
    if v >= 0:
        if flag_max and not flag_min:
            lamb = 1.0/(1.0 + np.exp(-r_pos*(v - v_off_pos)))
        else:
            lamb = 1.0/(1.0 + np.exp(-r_pos*(v - v_threshold_pos)))
        if v == MAX:
            flag_max = True
            flag_min = False
    elif v < 0:
        if flag_max and not flag_min:
            lamb = 1.0/(1.0 + np.exp(r_neg*(v - v_threshold_neg)))
        else:
            lamb = 1.0/(1.0 + np.exp(r_neg*(v - v_off_neg)))
        if v == MIN:
            flag_max = False
            flag_min = True
    i0 = i0_min + (i0_max - i0_min)*lamb*.01
    # Rs and alpha are constants in this approach
    inside = alpha*i0*Rs*d*np.exp(alpha*d*(abs(v) + Rs*i0))
    counter = 0
    w = 110
    aux = w - 1 # not equal to w MUST
    while round(aux,5) != round(w,5):  # tipico procedicimiento para resolwer lambert
#    while counter < 100:
        aux = w
        # Newton-Raphson
        w = w + (inside*np.exp(-w) - w)/(1.0 + w)
        counter += 1
#        print "Iteration", counter, "w", w
    w=aux
    current =np.absolute( np.sign(v)*(w/(alpha*Rs) + d*(G1*abs(v)-i0) + G2*abs(v)))
    current1 =np.absolute( np.sign(v)*(w/(alpha*Rs) + d*(G1*abs(v)-i0) + G2*abs(v)))+abs(v)/2e13+5e-13
    # storing in arrays    
    lamb_array.append(lamb)
    i0_array.append(i0)
    w_array.append(w)
    current_array.append(current)
    current_array1.append(current1)
    voltage_array.append(v)

# plot all data
graph = plt.figure(1)
#plt.semilogy(voltage_array, current_array,  'bo-', label = 'original')
plt.semilogy(voltage_array, current_array1, '-', label = 'ohmic_leakage')
plt.semilogy(voltage_data, current_data, 'go-.', label = 'experimental')
plt.legend(loc = 'best')
plt.grid(True)
graph = plt.figure(2)
plt.plot(voltage_array, i0_array, 'ro-', label = 'i0')
plt.legend(loc = 'best')
plt.grid(True)
graph = plt.figure(3)
plt.plot(voltage_array, w_array, 'go-', label = 'lambert')
plt.legend(loc = 'best')
plt.grid(True)
graph = plt.figure(4)
plt.plot(voltage_array, lamb_array, 'ko-', label = 'lambda')
plt.legend(loc = 'best')
plt.grid(True)

f = open ("fichero0553.dat", "w")


for i in range(len(voltage_array)):
    f.write("%5.18f %5.18f\n" % (voltage_array[i], current_array[i]))


f.close()

print ('End of simulation')