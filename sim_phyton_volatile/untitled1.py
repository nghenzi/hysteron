# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 14:31:23 2019

@author: titi
"""

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
v_off_pos=3
v_threshold_neg = -7.5
v_off_neg = -3

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


Lmax=1
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
#            lamb = Lmax*1.0/(1.0 + np.exp(-r_pos*(v - v_off_pos/beta_th)))
            lamb = Lmax*np.tanh(r_pos*(v - v_off_pos/beta_th))/2+.5
        else:
            lamb = Lmax*np.tanh(r_pos*(v - v_threshold_pos))/2+.5
        if v == MAX:
            flag_max = True
            flag_min = False
    elif v < 0:
        alpha=alpha_neg
        i0_max=i0_max_neg*Lmax
        if flag_max and not flag_min:
            lamb = Lmax*np.tanh(-r_neg*(v - v_threshold_neg))/2+.5
        else:
            lamb = Lmax*np.tanh(-r_neg*(v - v_off_neg/beta_th))/2+.5
        if v == MIN:
            flag_max = False
            flag_min = True
    lamb_array.append(lamb)
    
    
    
plt.plot(v_sweep, lamb_array, 'r-', linewidth=3.0,  label = 'simulation')