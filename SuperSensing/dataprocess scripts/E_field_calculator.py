# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 15:13:22 2023

@author: Admin
"""
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot

fname = 'C:/data/2023/2-16-2023/export_F08_full_reflection_1.dat'

spec = np.loadtxt(fname)

t = spec[:,0]
dI = spec[:,1]
c = 3e8
r_41 = 3.9*1e-12
n_ir = 2.78
I_p = 0.45e-3/pi/(1.5e-3)**2
w = 2*pi*c/1030e-9
L = 1e-3
a = 80/180*pi
phi = a
E_thz = dI*2*c/I_p/(w*n_ir**3*r_41*L)/(np.cos(a)*np.sin(2*phi)+2*np.sin(a)*np.cos(2*phi))

const1 = 2*c/I_p/(w*n_ir**3*r_41*L)
# print(const1)

figure, ax = plt.subplots(2,1)
ax1 = ax[0]
ax2 = ax[1]
ax1.plot(t, E_thz, label = 'estimated field stregnth')
ax1.set_xlabel('time (ps)')
ax1.set_ylabel('E (V/m)')
ax1.legend(loc = 'best')

ax2.plot(t, dI, label = 'original signal')
ax2.set_xlabel('time (ps)')
ax2.set_ylabel('I (V)')
ax2.legend(loc = 'best')