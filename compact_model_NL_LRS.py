# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 15:11:15 2021

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt

i0_max, i0_min = 220e-6, 1e-8
vth_pos, vth_neg = 2.8, 1.4
r_pos, r_neg = 0.5, 0.3

e = float('inf')

alpha_pos, alpha_neg = 2.1, 1.2
alpha, Rs = 3.9, 0.2e3


MAX, MIN = 3.0, -3.0
Npoints = 500 

x = np.arange(0,3,0.0001)
v = np.concatenate((x,x[::-1]))#,-x,-x[::-1]))

fig, ax = plt.subplots()

for idx in range(1,20):
    # R2 = 10**(R2)
    alpha = 0.1+0.2*idx
    
    Rs = 1
    i0_max = 100e-6
    R2 = float('inf')
    Rs*=1.
    R1, R2 = 1e7, R2
    G1, G2 = 1.0/R1, 1.0/R2
    d = 1.0/(1 + Rs*G1)

    l = np.concatenate((np.tanh((x-vth_pos)/r_pos)/2+.5,[1]*x.shape[0]))#,
                        # np.tanh((-x+vth_neg)/r_neg)/2+.5, [0]*x.shape[0]))
    i0 = i0_min + (i0_max - i0_min)*l
    jk = alpha*i0*Rs*d*np.exp(alpha*d*(abs(v) + Rs*i0))
    w = np.log(1+jk)*( 1- ( (np.log(1+(np.log(1+jk))))/(2+(np.log(1+jk))) ) )
            
    current = np.absolute( np.sign(v)*(w/(alpha*Rs) + d*(G1*abs(v)-i0) + G2*abs(v)))
    
    print (idx,alpha,  plt.cm.RdYlGn(idx/15) )
    ax.plot(v[3001:],current[3001:], linewidth= 2, color = 'k' )
# ax.set_xlim(-5,5)
# ax.set_ylim(1e-7,0.3e-3)

# ax.set_xlabel('V [V]', fontsize=16)
# ax.set_ylabel('I [A]', fontsize=16)

# ax.tick_params(direction='in')

# from matplotlib.ticker import MultipleLocator
# ax.xaxis.set_minor_locator(MultipleLocator(0.4))
# ax.xaxis.set_tick_params(which='minor', direction='in')
# ax.yaxis.set_tick_params(which='minor', direction='in')
                                           
# plt.setp(ax.get_xticklabels(), fontsize=14)
# plt.setp(ax.get_yticklabels(), fontsize=14)

# plt.tight_layout()
#%%
for i in range(10) : ax.lines[-1].remove()
# for line in ax.lines:
#     line.set_color('grey')
#%%
ax.set_yscale('linear')
ax.set_yscale('log')
ax.set_ylim(1e-7,5e-3)
