# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 09:08:40 2023

@author: Admin
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 08:32:29 2023

@author: Admin
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 08:54:05 2023

@author: Admin
"""
from numba import njit
import numpy as np
import dataprocessor as dp
import matplotlib.pyplot as plt
import plottool_v6 as pt
from numpy import pi
from scipy.signal import stft, spectrogram
import glob
import re
import tkinter as tk
import Fast_data_process as fdp
from numpy import pi


def fftx_blackman(t, x, pad = 2):
    ts = abs(t[1]-t[0])
    NFFT = int(2**(np.ceil(np.log2(len(t)))+pad))
    fs = 1/ts
    x_reduce = x-np.average(x)
    x_filter = x_reduce*np.blackman(len(x_reduce))
    SX = np.fft.fft(x_filter,n=NFFT)
    sx = SX[0:int(NFFT/2)]
    freq = fs/2*np.linspace(0,1,len(sx))
    return freq, sx

def windowed_energy_trans(idx1, idx2, x_sam, x_ref):
    I_sam = np.trapz(x_sam[idx1:idx2]**2)/np.trapz(x_ref[idx1:idx2]**2)
    return I_sam

path = 'C:/Users/Admin/Dropbox/Data/SPEC SO/2023/Aug/08-03-2023'

energy_trans = np.zeros(13)
energy_trans_so = np.zeros(13)
delays = np.arange(0, 1.21, step = 0.1)/1.2*2

for i in range(0,13):
    # plt.figure()
    # plt.xlabel('time (ps)')
    # plt.ylabel('Amplitude (a.u.)')
    t_ref_raw, x_ref_raw, y_ref_raw = fdp.load_single_spec(path, 'Delay{}_D0_IFE'.format(i))
    t_sam_raw, x_sam_raw, y_sam_raw = fdp.load_single_spec(path, 'Delay{}_D1_IFE'.format(i))
    t_sam = t_sam_raw[:, 0]
    idx0 = int(len(t_sam)/2)
    idx1 = np.where(t_sam>=37)[0][0]
    idx2 = np.where(t_sam<=38.4)[0][-1]
    x_sam = np.average(x_sam_raw, axis = 1)
    # x_sam = x_sam - np.mean(x_sam)
    x_ref = np.average(x_ref_raw, axis = 1)
    # x_ref = x_ref - np.mean(x_ref)
    # plt.plot(t_sam, x_sam, label = 'evanescent')
    # plt.plot(t_sam, x_ref, label = 'reference')
    energy_trans[i] = windowed_energy_trans(0, -1, x_sam, x_ref)
    energy_trans_so[i] = windowed_energy_trans(idx1, idx2, x_sam, x_ref)   

# plt.figure()
# # plt.plot(t_so_ref, x_so_ref, label = 'reference')
# plt.xlabel('time (ps)')
# plt.ylabel('Amplitude (a.u.)')
# t_sam, x_sam, y_sam = fdp.load_single_spec(path, 'D5_SO')
# t_so_sam = t_sam[:, 0]
# x_so_sam = np.average(x_sam, axis = 1)
# x_so_d4 = x_so_sam - np.mean(x_so_sam)
# t_ref, x_ref, y_ref = fdp.load_single_spec(path, 'D0_SO')
# t_so_ref = t_ref[:, 0]
# x_so_ref = np.average(x_ref, axis = 1)
# x_so_ref = x_so_ref - np.mean(x_so_ref)
# idx0 = int(len(t_so_sam)/2)


# energy_trans['SO'] = x_so_sam**2/x_so_ref**2

# plt.plot(t_so_sam, x_so_sam, label = 'evanescent')
# plt.plot(t_so_sam, x_so_ref, label = 'reference')
# plt.legend()
# plt.grid(1)


# time_plot = t_so_ref




plt.figure()

# plt.plot(delays, energy_trans, label = 'complete energy transmit', marker = 'x')
plt.plot(delays, energy_trans_so, label = 'IFE region energy transmit', marker = 'x')
plt.xlabel('phase delays ($\pi$)')
plt.ylabel('Transmitted Intensity (a.u.)')
# plt.xlim((0,1.75))
# plt.ylim((0,0.05))
# plt.legend()
# plt.grid(1)
# for key in list(energy_trans):
#     plt.plot(delays, energy_trans[key], label = key)
#     plt.xlabel('center of time window (ps)')
#     plt.ylabel('Energy Transmitted (a.u.)')
#     plt.grid(1)





    

    
    
    
    
    