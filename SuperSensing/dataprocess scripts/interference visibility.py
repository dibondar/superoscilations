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

path = 'C:/Users/Admin/Dropbox/Data/SPEC SO/2023/Jul/07-06-2023'
all_fnames = fdp.get_fnames(path)
x = {}
y = {}
t = {}
zxx = {}
f = {}
t_stft = {}
fnames_delay = []

f_ave_all = {}
for file in all_fnames:
    if 'IFE' in file and 'Delay' in file:
        fnames_delay.append(file)
fnames_delay.sort(key = fdp.str_get_Delay)
for fname in fnames_delay:
    t[fname], x[fname], y[fname] = fdp.load_single_spec(path, fname)
t_ave, x_ave = fdp.average_spec_noref(t, x)

t_ref, x_ref, y_ref = fdp.load_single_spec(path, 'D0_F08')
t_f08 = t_ref[:, 0]
x_f08 = np.average(x_ref, axis = 1)
t_ref, x_ref, y_ref = fdp.load_single_spec(path, 'D0_F07')
t_f07 = t_ref[:, 0]
x_f07 = np.average(x_ref, axis = 1)
# f0, sx0 = fftx_blackman(t0, x0)

# intensity = np.zeros(len(list(x_ave)))
inter_vis = np.zeros(len(list(x_ave)))
energy_trans = np.zeros(len(list(x_ave)))
for n,key in enumerate(list(x_ave)):
    twindow = t_ave[key][(t_ave[key]>32)&(t_ave[key]<36)]
    xwindow = x_ave[key][(t_ave[key]>32)&(t_ave[key]<36)]
    x_f08_shift = fdp.spec_shift(t_ave[key], x_ave[key], 0.1*(n+1))
    x0_ife = (x_f07+x_f08_shift)[(t_ave[key]>32)&(t_ave[key]<36)]
    # plt.figure()
    # f, sx = fftx_blackman(t_ave[key], x_ave[key])
    # transmission = abs(sx)/abs(sx0)
    energy_trans[n] = np.trapz(xwindow**2, twindow)/np.trapz(x0_ife**2, twindow)
    # trans_pick = sx[(f>0.6)&(f<0.9)]
    # inter_vis[n] = fdp.interferometric_vis(trans_pick)
    # plt.plot(f, abs(sx))
    # plt.xlim((0, 2))
    # f_so = f_ave[(t_stft[key] >= 4.4)&(t_stft[key] <= 5.6)]
    # f_so_ave[n] = np.trapz(f_so, dx = 0.05)/1.2
delay = np.arange(0.1, 1.21, step = 0.1)

plt.figure()
plt.plot(delay, energy_trans, label = 'experiment')
plt.xlabel('applied time delay (ps)')
plt.ylabel('Energy Transmitted (a.u.)')


# root = tk.Tk()
# plotgo = pt.Plottool(root, delay, inter_vis)
# root.mainloop()

time = np.arange(0, 10.01, step = 0.05 )
# x1 = np.sin(2*pi*0.8*time)
x2 = fdp.generate_gauss_cos(time, 0.7, sigma = 5)
x1 = fdp.generate_gauss_cos(time, 0.8, sigma = 5)
f0_theo, sx0_theo = fftx_blackman(time, (x1+x2))


inter_vis_theo = np.zeros(len(delay))
energy_trans_theo = np.zeros(len(delay))
for i in range(0, 12):
    time_delay = 0.1*(i+1)+time
    x_delay = fdp.generate_gauss_cos(time_delay, 0.8, sigma = 5)
    x_ife = (x2 + x_delay)[(time>2)&(time<6)]
    time_window = time[(time>2)&(time<6)]
    f, sx = fftx_blackman(time, x_ife)
    transmission = abs(sx)/abs(sx0_theo)
    energy_trans_theo[i] = np.trapz(x_ife**2, time_window)/np.trapz(((4*x1+4*x2)**2)[(time>2)&(time<6)], time_window)
    # plt.figure()
    # plt.plot(time, x_ife)
    # plt.xlim((0, 10))
    trans_pick = transmission[(f>0.6)&(f<0.9)]
    inter_vis_theo[i] = fdp.interferometric_vis(trans_pick)

    
# plt.figure()
# plt.plot(delay, energy_trans_theo, color = 'red', label = 'calculate')
plt.legend()
# plt.xlabel('time delay (ps)')
# plt.ylabel('Interferometric visibility (a.u.)')

# F07 = path+'D1_F07.dat'
# E_f07 = path
# plt.figure()
# plt.plot(time, )


    

    
    
    
    
    