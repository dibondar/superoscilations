
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 14:34:48 2023

@author: Admin
"""

import numpy as np
import dataprocessor as dp
import matplotlib.pyplot as plt
from numpy import pi
import glob
from scipy.integrate import simps
from scipy.signal import stft, spectrogram
from numba import njit
# pad = np.zeros(100)
# t = np.arange(0, 1.2, step = 0.05)
# x = np.sin(2*pi*2*t)
# f, sx = dp.Basic_functions().array_fftx(x, t)
# fig, ax = plt.subplots(1,1,figsize = (9,6), dpi=300)
# ax.plot(f, abs(sx))

def load_single_spec(folder, fname, ftype = '.dat'):
    folder = folder.replace('\\','/')
    flist = glob.glob(folder+'/'+fname+'_*[0-9]'+ftype)
    t_store = []
    x_store = []
    y_store = []
    for i in range(0, len(flist)):
        fnum = str(i+1)
        try:
            data = np.loadtxt(folder+'/'+fname+'_'+fnum+ftype)
        except OSError:
            continue
        else:
            data = np.loadtxt(folder+'/'+fname+'_'+fnum+ftype)
            t_store.append(data[:,0])
            x_store.append(data[:,1])
            y_store.append(data[:,2])
    t_store = np.array(t_store).T
    x_store = np.array(x_store).T
    y_store = np.array(y_store).T
    return t_store, x_store, y_store


def fftx_stft(t, x, pad = 2):
    ts = abs(t[1]-t[0])
    NFFT = int(2**(np.ceil(np.log2(len(t)))+pad))
    fs = 1/ts
    x_normal = x-np.average(x)
    # freq_auto = 0
    # t_auto = 0
    # zx = 0
    freq_auto, t_auto, zx = stft(x_normal, fs, nperseg = 100, noverlap=99,
                                        window = ('gaussian',14), nfft = NFFT)
    return freq_auto, t_auto, zx

path = 'C:/Users/Admin/Dropbox/Data/SPEC SO/2023/Jun/06-30-2023'
t_all, x_all, y_all = dp.Basic_functions.load_single_spec(path, 'D0_SO')

t_all, x_all, y_all = load_single_spec(path, 'D0_SO')

time = t_all[:, 0]
time_mid = round((time[0]+time[-1])/2)
time_diff = time_mid - 0
time = time - time_diff
x = np.average(x_all, axis = 1)
f, t, Zxx = fftx_stft(time, x)
nrow, ncol = np.shape(Zxx)
f_ave = []
for i in range(0, ncol):
    ampl = abs(Zxx)
    ampl_sum = sum(ampl[:,i])
    f_ave.append(sum(f*ampl[:,i])/ampl_sum)
    
f_ave = np.array(f_ave)
# f_so_ave = np.mean(f_ave[(time >= -0.6) & (time <= 0.6)])
# plt.plot(time[(time >= -0.6) & (time <= 0.6)], f_ave[(time >= -0.6) & (time <= 0.6)])
plt.plot(t, f_ave)
plt.xlabel('time (ps)')
plt.ylabel('Instant frequency(THz)') 







