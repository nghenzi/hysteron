import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
import pandas as pd 

plt.close('all')
#
filename="volatile.dat"
data_100nA = np.loadtxt(filename, delimiter="\t", skiprows=2, usecols=[0,1,2])
data_10nA = np.loadtxt(filename, delimiter="\t", skiprows=2, usecols=[6])