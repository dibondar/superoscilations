
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 09:58:25 2023

@author: Admin
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 10:55:04 2022

@author: ppsapphire
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy.signal import hilbert
import dataprocessor as dp




bf = dp.Basic_functions()
# path = '/home/ppsapphire/Dropbox/SO project/artificial pulses'
path = 'C:/Users/Admin/Dropbox/SO project/artificial pulses'
# path = '/Users/ppsap/Dropbox/SO project/artificial pulses'


'''read in pulses'''

filename1 = 'Glucose_calculated_1.dat'
filename2 = 'LGA_calculated_1.dat'
spec1 = np.loadtxt(path + '/' + filename1)
spec2 = np.loadtxt(path + '/' + filename2)
# spec2 = np.loadtxt(path + '/' + filename2)
t0 = 15.00
t1 = np.around(spec1[:,0], decimals = 2)
t_shift = t1[0] - t0
t2 = spec2[:,0]
x1 = spec1[:,1]
t_new = t1 - t_shift
x2 = spec2[:,1]
zerofiller = np.zeros_like(t1)

new_spec1 = np.stack((t_new, x1, zerofiller), axis = 1)
new_spec2 = np.stack((t_new, x2, zerofiller), axis = 1)

new_fname1 = 'Glucose_cal_shifted_1.dat'
new_fname2 = 'LGA_cal_shifted_1.dat'

np.savetxt(path + '/' + new_fname1, new_spec1, fmt = '%1.8f')
np.savetxt(path + '/' + new_fname2, new_spec2, fmt = '%1.8f')
# x2 = spec2[:,1]
plt.figure()
plt.plot(t_new, x1, t_new, x2)








