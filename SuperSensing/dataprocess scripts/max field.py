# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 08:54:05 2023

@author: Admin
"""
import numpy as np
import dataprocessor as dp
import matplotlib.pyplot as plt
import plottool_v6 as pt
from numpy import pi
from scipy.signal import stft, spectrogram
import glob
import re
import tkinter as tk

def get_fnames(path, filetype='.dat'):
    pathlist = glob.glob(path+'/'+'*'+filetype)
    filenames = []
    for i in range(0,len(pathlist)):
        pathlist[i] = pathlist[i].replace('\\','/')
        iter1 = re.finditer('/', pathlist[i])
        iter2 = re.finditer('_', pathlist[i])
        for j in iter1:  #loop to last / sign
            id1 = j
        for k in iter2: #loop ot last _ sign
            id2 = k
        id_backs = id1.end() #get the index after last /
        id_unders = id2.start() # get the index before last _
        name = pathlist[i][id_backs:id_unders]
        if name in filenames:
            continue
        else:
            filenames.append(name)
    return filenames



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



def average_spec(t, x):
    x_ave = {}
    t_ave = {}
    for key in list(x):
        x_ave[key] = np.average(x[key], axis = 1)
        t_ave[key] = t[key][:, 0]
    return t_ave, x_ave
        
def str_get_D(string):
    str_seg_all = string.split('_')
    for str_seg in str_seg_all:
        if 'D' in str_seg:
            pick_str = str_seg
            continue
    num = ''
    for char in pick_str:
        if char.isdigit():
            num = num + char
    num = int(num)
    return num

def fftx_stft(t, x, pad = 2):
    ts = abs(t[1]-t[0])
    NFFT = int(2**(np.ceil(np.log2(len(t)))+pad))
    fs = 1/ts
    x_normal = x-np.average(x)
    phi0 = np.arcsin(x[0]/max(x_normal))/2/pi/0.8
    phi1 = np.arcsin(x[1]/max(x_normal))/2/pi/0.8
    art_t0 = np.arange(phi0-5, phi0, step = 0.05)
    art_t1 = np.arange(phi1, phi1+5, step = 0.05)
    pad0 = max(x_normal)*np.sin(0.8*2*pi*art_t0)
    pad1 = max(x_normal)*np.sin(0.8*2*pi*art_t1)

    x_pad = np.concatenate((pad0, x_normal, pad1))
    freq_auto, t_auto, zx = stft(x_pad, fs, nperseg = len(t)/2, noverlap=len(t)/2-1,
                                        window = ('gaussian',12), nfft = NFFT)
    return freq_auto, t_auto, zx

def interferometric_vis(x):
    i_max = max(x**2)
    i_min = min(x**2)
    ivis = (i_max-i_min)/(i_max+i_min)
    return ivis

path = 'C:/Users/Admin/Dropbox/Data/SPEC SO/2023/Jun/06-30-2023'
all_fnames = get_fnames(path)
x = {}
y = {}
t = {}
zxx = {}
f = {}
t_stft = {}
fnames_so = []

f_ave_all = {}
for file in all_fnames:
    if 'SO' in file:
        fnames_so.append(file)
fnames_so.sort(key = str_get_D)
for fname in fnames_so:
    t[fname], x[fname], y[fname] = load_single_spec(path, fname)
t_ave, x_ave = average_spec(t, x)
x_max = []
f_so_ave = np.zeros(len(list(x_ave)))
intensity = np.zeros(len(list(x_ave)))
inter_vis = np.zeros(len(list(x_ave)))
plt.figure()
for n,key in enumerate(list(x_ave)):
    x_max.append(np.max(x[key]))
    f[key], t_stft[key], zxx[key] = fftx_stft(t[key][:,0], x_ave[key])
    nrow, ncol = np.shape(zxx[key])
    ampl = abs(zxx[key])
    f_ave = np.zeros(ncol)
    x_so = x_ave[key][(t[key][:,0] >= 34.4)&(t[key][:,0] <= 35.6)]
    for i in range(0, ncol):
        ampl_sum = sum(ampl[:,i])
        f_ave[i] = (sum(f[key]*ampl[:,i])/ampl_sum) 
    f_ave_all[key] = f_ave
    plt.plot(t_stft[key], f_ave, label = key)
    
    intensity[n] = np.trapz(x_ave[key]**2, dx = 0.05)
    # inter_vis[n] = interferometric_vis(x_so)
    # f_so = f_ave[(t_stft[key] >= 4.4)&(t_stft[key] <= 5.6)]
    # f_so_ave[n] = np.trapz(f_so, dx = 0.05)/1.2
plt.legend('on')
distance = np.array([0,13.5, 27, 99, 159, 249, 339, 489])
# distance = np.array([0,13.5, 27, 99, 159, 249, 339, 489, 1989, 3489, 4989])
# intensity_FHs = {}
plt.figure()
for m in range(0, 4):
    fnames_fhs = []
    for file in all_fnames:
        if 'F0{}'.format(str(m+5)) in file:
            fnames_fhs.append(file)
    fnames_fhs.sort(key = str_get_D)
    for fname in fnames_fhs:
        t[fname], x[fname], y[fname] = load_single_spec(path, fname)
    t_ave, x_ave = average_spec(t, x)
    x_max = []
    # f_so_ave = np.zeros(len(list(x_ave)))
    intensity = np.zeros(len(fnames_fhs))
    for n,key in enumerate(fnames_fhs):
    #     x_max.append(np.max(x[key]))
    #     f[key], t_stft[key], zxx[key] = fftx_stft(t[key][:,0], x_ave[key])
    #     nrow, ncol = np.shape(zxx[key])
    #     ampl = abs(zxx[key])
    #     f_ave = np.zeros(ncol)
        # x_so = x_ave[key][(t[key][:,0] >= 34.4)&(t[key][:,0] <= 35.6)]
        # for i in range(0, ncol):
        #     ampl_sum = sum(ampl[:,i])
        #     f_ave[i] = (sum(f[key]*ampl[:,i])/ampl_sum) 
        # f_ave_all[key] = f_ave
        intensity[n] = max(x_ave[key]**2)
    plt.semilogy  (distance/300, intensity, label = 'intensity_F0{}'.format(str(m+5)))

plt.xlabel('d ($\lambda$)')
plt.ylabel('max intensity')
plt.legend()
    # intensity_FHs['intensity_F0{}'.format(str(m+5))] = intensity



    

    
# x_max = np.array(x_max)

# distance = np.array([0,13.5, 27, 99, 159, 249, 339, 489])
# root = tk.Tk()
# plotgo = pt.Plottool(root, t_stft, f_ave_all)
# root.mainloop()


# plt.figure()
# ax = plt.axes()
# ax1 = ax.twinx()
# ax.semilogy(distance / 300, inter_vis, '*-')
# #ax.semilogy(distance / 300, intensity, color = 'red')
# # ax1.semilogy(distance / 300, intensity, '*-', color = 'red')
# ax.set_xlabel('$\lambda$ (um)')
# ax.set_ylabel('Overall max amplitude (V)')
# ax1.set_ylabel('Intensity', color = 'red')
# ax.set_xlim((0, 1))


    

    
    
    
    
    