# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 11:21:26 2023

@author: Admin
"""

import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import copy
from statistics import mean
from scipy.signal.windows import blackman
from numpy import pi
from scipy.signal import hilbert
from numba import njit
from numba.experimental import jitclass

@njit    
def formatinput(x):
    if type(x) is list:
        x_new = dict()
        for i in range(len(x)):
            name = 'x'+str(i)
            x_new[name] = x[i]
    if type(x) is np.ndarray:
        x_new = dict()
        if np.size(x) == len(x):
                name = 'x0'
                x_new[name] = x
        else:
            for i in range(np.size(x,1)):
                name = 'x'+str(i)
                x_new[name] = x[:,i]
    if type(x) is dict:
        x_new = x
    return x_new

def array_chop_pad(x, y, x0, x1):
    idx0 = np.where(x >= x0)[0][0]
    idx1 = np.where(x <= x1)[0][-1]
    pad0 = np.zeros_like(x[0:idx0]) + y[idx0]
    pad1 = np.zeros_like(x[idx1:]) + y[idx1]
    y_new = np.concatenate((pad0, y[idx0:idx1], pad1), axis = 0)
    return y_new

def array_chop(x, y, x0, x1):      
    # x0 = round(x0, ndigits = 2)
    # x1 = round(x1, ndigits = 2)
    # step = round(abs(x[1]-x[0]), ndigits = 2)
    idx0 = np.where(x >= x0)[0][0]
    idx1 = np.where(x <= x1)[0][-1]
    x_new = x[idx0:idx1+1]
    y_new = y[idx0:idx1+1]
    return x_new, y_new

def findzeropoints(x,y):
    for i in range(0,len(x)-1):
        x_zeros = []
        if y[i]*y[i+1] <= 0:
           x_zeros.append((x[i]+x[i+1])/2)
    return np.array(x_zeros)

def find_1st_zero(x,y):
    x0 = None
    for i in range(0,len(y)-1):
        if y[i]*y[i+1] <= 0:
           x0 = (x[i]+x[i+1])/2
        if i == len(y)-2 and x0 == None:
            x0 = x[round(len(x)/2)]
    return x0

def normalize_signal(x):
    x = formatinput(x)
    x_normal = {}
    keys = list(x)
    for key in keys:
        x_normal[key] = np.zeros_like(x[key])
        dim = np.shape(x[key])
        if len(dim) == 1:
            x_max = x[key].max()
            x_normal[key] = x[key]/x_max
        if len(dim) == 2:
            x_max = x[key].max(axis = 1)
            for i in range(0, dim[1]):
                x_normal[key][:, i] = x[key][:, i]/x_max[i]
    return x_normal

def normalize_signal_samref(x):
    x = formatinput(x)
    x_normal = {}
    keys = list(x)
    for key in keys:
        if key+'_ref' in keys:
            x_max = x[key+'_ref'].max()
            x_normal[key] = x[key]/x_max
        else:
            x_max = x[key].max()
            x_normal[key] = x[key]/x_max
    return x_normal

# def signal_noise_filter( x):
#     x = formatinput(x)
#     x_filter = {}
#     keys= list(x)
#     for key in keys:
#         x_filter[key] = np.zeros_like(x[key])
#         dim = np.shape(x[key])
#         if len(dim) == 1:
#             x_max = x[key].max()
#             for j in range(0, len(x[key])):
#                 x_normal = x[key][j]/x_max
#                 if x_normal <= 1e-3:
#                     x_filter[key][j] = 0
#                 else:
#                     x_filter[key][j] = x[key][j]
#         if len(dim) == 2:
#             x_max = x[key].max(axis = 1)
#             for i in range(0, dim[1]):
#                 for j in range(0, len(x[key][:, i])):   
#                     x_normal = x[key][j, i]/x_max[i]
#                     if x_normal <= 1e-3:
#                         x_filter[key][j, i] = 0
#                     else:
#                         x_filter[key][j, i] = x[key][j, i]
#     return x_filter

# def getindex(freq,E_sam,E_ref,L,E_echo=None,c=3e8):
#     freq = formatinput(freq)
#     E_sam = formatinput(E_sam)
#     E_ref = formatinput(E_ref)
#     fname = list(freq)[0]
#     samname = list(E_sam)[0]
#     refname = list(E_ref)[0]
#     omega = freq[fname]*2*np.pi*1e12
#     if E_echo is None:
#         phi = np.unwrap((np.angle(E_sam[samname]/E_ref[refname])))
        
#         n = dict()
#         n[samname] = phi/omega/L*c+1
#     else:
#         Lr = np.linspace(0.9*L,1.1*L,1000)
#         dn = np.zeros(len(Lr))
#         E_echo = formatinput(E_echo)
#         echoname = list(E_echo)[0]
#         # phi_sam = np.unwrap((np.angle(E_sam[samname]/E_ref[refname])),discont=pi/4)
#         # phi_echo = np.unwrap((np.angle(E_echo[echoname]/E_ref[refname])),discont=pi/4)
#         phi_sam = (np.angle(E_sam[samname]/E_ref[refname]))
#         phi_echo = (np.angle(E_echo[echoname]/E_ref[refname]))
#         for i in range(0,len(Lr)):
#             n_sam = phi_sam/omega/Lr[i]*c+1
#             n_echo = (phi_echo/omega/Lr[i]*c+1)/3          
#             dn[i] = np.nanmean(n_echo-n_sam)
#         L0 = find_1st_zero(Lr,dn)
#         n = dict()
#         n[samname] = phi_sam/omega/L0*c+1
#     return n

def getindex_array(freq,sx_sam,sx_ref,L,SX_echo=None,c=3e8, ph_mod = None):
    '''
    Parameters
    ----------
    freq : np.array
        DESCRIPTION.frequency of the spectrums
    E_sam : 1d np.array
        DESCRIPTION.sample spectrum
    E_ref : 1d np.array
        DESCRIPTION.reference spectrum
    L : float
        DESCRIPTION.initial guess of the sample thickness
    E_echo : np.array, optional
        DESCRIPTION. echo spectrum
    c : float, optional
        DESCRIPTION.speed of light, The default is 3e8.
    Returns
    -------
    n : np.array
        DESCRIPTION. calcullated refractive index
    '''
    
    
    
    # E_sam = sx_sam
    # E_ref = sx_ref
    # freq, E_sam, E_ref = badtrans_remove(freq, sx_sam, sx_ref)
    freq, E_sam, E_ref = FD_noise_remove(freq, sx_sam, sx_ref)
    omega = freq*2*pi*1e12
    
    freq  = freq
    if SX_echo is None: #if there is no echo in the scan
        
        # phi = np.angle(E_sam/E_ref) #calculate phase of transmission
        if ph_mod == None:
            phi_sam = abs(np.unwrap(np.angle(E_sam), discont = pi/4))
            phi_ref = abs(np.unwrap(np.angle(E_ref), discont = pi/4))

        else:
            phi_sam = abs(np.unwrap(np.angle(E_sam), discont = pi/4))
            phi_ref = abs(np.unwrap(np.angle(E_ref), discont = pi/4))
            phi_sam = phi_sam + ph_mod
            # phi_ref += ph_mod
        
        phi = np.unwrap(phi_sam - phi_ref, discont = pi/4)
        # phi = np.unwrap(phi, discont = pi/4)
        # phi = abs(phi)
        n = phi/omega/L*c+1 # calculate refractive index
        # n = index_normal(freq, n)
        
    else:
        freq, E_sam, E_ref, E_echo = badtrans_remove(freq, sx_sam, sx_ref, E_echo = SX_echo)
        Lr = np.linspace(0.8*L,1.2*L,1000) #create a range of estimated thickness
        dn = np.zeros((len(freq),len(Lr))) #prepare dn to store the refractive index difference
        # L0 = np.zeros(len(freq))
        # dn_max = np.zeros(len(Lr))
        # dn_min = np.zeros(len(Lr))
        phi_sam = np.unwrap((np.angle(E_sam/E_ref)),discont=pi/4) # phase angle of transmission from sample signal
        phi_echo = np.unwrap((np.angle(E_echo/E_ref)),discont=pi/4)#phase angle of transmission from echo signal
        # phi_sam = (np.angle(E_sam/E_ref))
        # phi_echo = (np.angle(E_echo/E_ref))
        for i in range(0,len(Lr)):
            n_sam = phi_sam/omega/Lr[i]*c+1  #refractive index calculated from sample signal
            # n_sam = index_normal(n_sam)
            n_echo = (phi_echo/omega/Lr[i]*c+1)/3 #refractive index calculated from echo signal
            # n_echo = index_normal(n_echo)
            dn[:,i] = n_echo-n_sam # difference of refractive index average over entire frequency range
            # dn_max[i] = max(n_echo-n_sam)
            # dn_min[i] = min(n_echo-n_sam)
        # L0 = find_1st_zero(Lr,dn) #get the true thickness by finding the zeros point
        L0 = L
        n = phi_sam/omega/L0*c+1 #use the true thickness to calculate refractive index
        # n = index_normal(n)
        
    a = -2/L*np.log(abs(E_sam)/abs(E_ref)*(1+n)**2/4/n)
    # k = np.imag(n)
    # a = -4*pi*freq*k/c
    return  freq, E_sam, E_ref, n, a, phi_sam, phi_ref

def badtrans_remove( freq, E_sam, E_ref, freq_limit = 3):
    idx_limit = np.where(freq <= freq_limit)[0][-1]
    T = (abs(E_sam)/abs(E_ref))[0:idx_limit]
    idxs = np.where(T >= 1)[0]
    d0 = 0
    idx0 = 0
    for i in range(0, len(idxs)-1):
        d1 = idxs[i+1]-idxs[i]
        if d1 >= d0:
            idx0 = idxs[i] + 1
            idx1 = idxs[i+1] - 1
            d0 = d1
    # idx0 = idx0 
    # idx1 = idx1
    E_sam = E_sam[idx0:idx1]
    E_ref = E_ref[idx0:idx1]
    freq = freq[idx0:idx1]
    return freq, E_sam, E_ref

def FD_noise_remove(f, sx_sam, sx_ref, sx_echo = None):
    idx0 = np.where(f>0.1)[0][0]
    idx1 = np.where(f>3)[0][0]
    # idx1 = np.where(f<8)[0][-1]
    noise_sam = np.max(abs(sx_sam[idx1:]))
    # noise_ref = np.max(abs(sx_sam[idx0:]))
    idx_sam_0 = idx0
    idx_sam_1 = np.where(abs(sx_sam)>noise_sam)[0][-1]
    sx_sam_filter = sx_sam[idx_sam_0:idx_sam_1]
    
    # idx_ref_0 = np.where(abs(sx_ref)>noise_ref)[0][0]
    # idx_ref_1 = np.where(abs(sx_ref)>noise_ref)[0][-1]
    sx_ref_filter = sx_ref[idx_sam_0:idx_sam_1]
    
    f_filter = f[idx_sam_0:idx_sam_1]
    return f_filter, sx_sam_filter, sx_ref_filter
    

def weak_signal_remove(x, threshold = 1e-2):
    # x = formatinput(x)
    for key in x:
        x_temp = x[key]
        if len(np.shape(x_temp)) == 1:
            max_x = max(x_temp)
            for i in range(0, len(x_temp)):
                ratio = abs(x_temp[i])/max_x
                if ratio <= threshold:
                    x_temp[i] = x_temp[i] / 10
        else:
            max_x = x_temp.max(axis = 1)
            for i in range(0, len(max_x)):
                for j in range(0, len(x_temp[:, i])):
                    ratio = abs(x_temp[j, i])/max_x[i]
                    if ratio <= threshold:
                        x_temp[j, i] = x_temp[j, i] / 10
        x[key] = x_temp   
    return x
# def index_normal( f, n):
#     idx = np.where(f>3)[0][0]
#     n_min = min(n[idx:])
#     if n_min < 1:
#         n = n+(1-n_min)
#     return n
# def unwrap_to_inf( phi):
#     for i in range(0,len(phi)-1):
#         if phi[i+1] - phi[i] >= pi:
#             phi[0:i+1] = phi[0:i+1] + 2*pi
#     return phi

# def reverse_unwrap( phi):
#     phi_r = np.flip(phi)
#     phi_wrap = np.unwrap(phi_r)
#     phi_out = np.flip(phi_wrap)
#     return phi_out

# def partial_unwrap( phi):
#     phi_out = np.zeros_like(phi)
#     phi_peak = max(phi)
#     idx = np.where(phi == phi_peak)[0][0]
#     phi_wrap = np.unwrap(phi[idx:])
#     phi_out[idx:] = phi_wrap
#     phi_out[0:idx] = phi[0:idx]
#     return phi_out
    


def fftx_filter(xvalues,tvalues,pad):
    xvalues = formatinput(xvalues)
    tvalues = formatinput(tvalues)
    sx = dict()
    freq = dict()
    for name in list(xvalues):
        dim = np.shape(xvalues[name])
        dim2 = np.shape(tvalues[name])
        if len(dim) == 1: 
            blackman_filter = blackman(len(xvalues[name]))
            if len(dim2) == 1:
                ts = abs(tvalues[name][1]-tvalues[name][0])
            if len(dim2) == 2:
                ts = abs(tvalues[name][1,0]-tvalues[name][0,0])
            NFFT = int(2**(np.ceil(np.log2(len(tvalues[name])))+pad))
            fs = 1/ts
            xvalue_reduce = xvalues[name]-np.average(xvalues[name])
            xvalue_filtered = xvalue_reduce*blackman_filter
            SX = np.fft.fft(xvalue_filtered,n=NFFT)
            sx[name] = SX[0:int(NFFT/2)]
            freq[name] = fs/2*np.linspace(0,1,len(sx[name]))
        else:
            blackman_filter = blackman(len(xvalues[name]))
            blackman_filter = blackman_filter[:,np.newaxis]
            if len(dim2) == 1:
                ts = abs(tvalues[name][1]-tvalues[name][0])
            if len(dim2) == 2:
                ts = abs(tvalues[name][1,0]-tvalues[name][0,0])
            NFFT = int(2**(np.ceil(np.log2(len(tvalues[name])))+pad))
            fs = 1/ts
            xvalue_reduce = xvalues[name]-np.average(xvalues[name],axis=0)
            xvalue_filtered = xvalue_reduce*blackman_filter
            SX = np.fft.fft(xvalue_filtered,n=NFFT,axis=0)
            sx[name] = SX[0:int(NFFT/2),:]
            freq[name] = fs/2*np.linspace(0,1,len(sx[name]))
            freq[name] = freq[name][np.where(freq[name]>0.1)]
            sx[name] = sx[name][np.where(freq[name]>0.1)]
    return freq,sx     


@njit
def fftx(xvalues,tvalues,pad):
    # xvalues = formatinput(xvalues)
    # tvalues = formatinput(tvalues)
    sx = {}
    freq = {}
    for name in list(xvalues):
        dim = np.shape(xvalues[name])
        dim2 = np.shape(tvalues[name])
        if len(dim) == 1: 
            if len(dim2) == 1:
                ts = abs(tvalues[name][1]-tvalues[name][0])
            if len(dim2) == 2:
                ts = abs(tvalues[name][1,0]-tvalues[name][0,0])
            NFFT = int(2**(np.ceil(np.log2(len(tvalues[name])))+pad))
            fs = 1/ts
            xvalue_reduce = xvalues[name]-np.average(xvalues[name],axis=0)
            SX = np.fft.fft(xvalue_reduce,n=NFFT)
            # sx[name] = SX
            sx[name] = SX[0:int(NFFT/2)]
            freq[name] = fs/2*np.linspace(0,1,len(sx[name]))
            # freq[name] = freq[name][np.where(freq[name]>0.1)]
            # sx[name] = sx[name][np.where(freq[name]>0.1)]
        else:
            if len(dim2) == 1:
                ts = abs(tvalues[name][1]-tvalues[name][0])
            if len(dim2) == 2:
                ts = abs(tvalues[name][1,0]-tvalues[name][0,0])
            NFFT = int(2**(np.ceil(np.log2(len(tvalues[name])))+pad))
            fs = 1/ts
            xvalue_reduce = xvalues[name]-np.average(xvalues[name],axis=0)
            SX = np.fft.fft(xvalue_reduce,n=NFFT,axis=0)
            sx[name] = SX[0:int(NFFT/2),:]
            # sx[name] = SX
            freq[name] = fs/2*np.linspace(0,1,len(sx[name]))
            # freq[name] = freq[name][np.where(freq[name]>0.1)]
            # sx[name] = sx[name][np.where(freq[name]>0.1)]
    return freq,sx 

def fftx_hilbert(xvalues,tvalues,pad):
    xvalues = formatinput(xvalues)
    tvalues = formatinput(tvalues)
    sx = dict()
    freq = dict()
    for name in list(xvalues):
        dim = np.shape(xvalues[name])
        dim2 = np.shape(tvalues[name])
        if len(dim) == 1: 
            if len(dim2) == 1:
                ts = abs(tvalues[name][1]-tvalues[name][0])
            if len(dim2) == 2:
                ts = abs(tvalues[name][1,0]-tvalues[name][0,0])
            NFFT = int(2**(np.ceil(np.log2(len(tvalues[name])))+pad))
            fs = 1/ts
            xvalue_reduce = xvalues[name] - np.average(xvalues[name],axis=0)
            xvalue_complex = hilbert(xvalue_reduce)
            SX = np.fft.fft(xvalue_complex,n = NFFT)
            # sx[name] = SX
            sx[name] = SX[0:int(NFFT/2)]
            freq[name] = fs/2*np.linspace(0,1,len(sx[name]))
            # freq[name] = freq[name][np.where(freq[name]>0.1)]
            # sx[name] = sx[name][np.where(freq[name]>0.1)]
        else:
            if len(dim2) == 1:
                ts = abs(tvalues[name][1]-tvalues[name][0])
            if len(dim2) == 2:
                ts = abs(tvalues[name][1,0]-tvalues[name][0,0])
            NFFT = int(2**(np.ceil(np.log2(len(tvalues[name])))+pad))
            fs = 1/ts
            xvalue_reduce = xvalues[name]-np.average(xvalues[name],axis=0)
            xvalue_complex = hilbert(xvalue_reduce)
            SX = np.fft.fft(xvalue_complex,n = NFFT, axis = 0)
            sx[name] = SX[0:int(NFFT/2),:]
            # sx[name] = SX
            freq[name] = fs/2*np.linspace(0,1,len(sx[name]))
            # freq[name] = freq[name][np.where(freq[name]>0.1)]
            # sx[name] = sx[name][np.where(freq[name]>0.1)]
    return freq,sx 

def fftx_omega(xvalues,tvalues,pad):
    xvalues = formatinput(xvalues)
    tvalues = formatinput(tvalues)
    sx = dict()
    omega = dict()
    for name in list(xvalues):
        dim = np.shape(xvalues[name])
        dim2 = np.shape(tvalues[name])
        if len(dim) == 1: 
            if len(dim2) == 1:
                ts = abs(tvalues[name][1]-tvalues[name][0])
            if len(dim2) == 2:
                ts = abs(tvalues[name][1,0]-tvalues[name][0,0])
            NFFT = int(2**(np.ceil(np.log2(len(tvalues[name])))+pad))
            fs = 1/ts
            SX = 1/2/np.pi*np.fft.fft(xvalues[name],n=NFFT)
            sx[name] = SX[0:int(NFFT/2)]
            omega[name] = 2*np.pi*fs/2*np.linspace(0,1,len(sx[name]))
        else:
            if len(dim2) == 1:
                ts = abs(tvalues[name][1]-tvalues[name][0])
            if len(dim2) == 2:
                ts = abs(tvalues[name][1,0]-tvalues[name][0,0])
            NFFT = int(2**(np.ceil(np.log2(len(tvalues[name])))+pad))
            fs = 1/ts
            SX = 1/2/np.pi*np.fft.fft(xvalues[name],n=NFFT,axis=0)
            sx[name] = SX[0:int(NFFT/2),:]
            omega[name] = 2*np.pi*fs/2*np.linspace(0,1,len(sx[name]))
    return omega,sx

@njit
def array_fftx(xvalues,tvalues,pad=2):
    ts = abs(tvalues[1]-tvalues[0])
    NFFT = int(2**(np.ceil(np.log2(len(tvalues)))+pad))
    fs = 1/ts
    xvalue_reduce = xvalues-np.average(xvalues)
    SX = np.fft.fft(xvalue_reduce,n=NFFT)
    sx = SX[0:int(NFFT/2)]
    freq = fs/2*np.linspace(0,1,len(sx))
    # freq = freq[np.where(freq>0.1)]
    # sx = sx[np.where(freq>0.1)]
    return freq,sx 

def array_fftx_hilbert(xvalues,tvalues,pad=2):
    ts = abs(tvalues[1]-tvalues[0])
    NFFT = int(2**(np.ceil(np.log2(len(tvalues)))+pad))
    fs = 1/ts
    xvalue_reduce = xvalues-np.average(xvalues)
    xvalue_hilbert = hilbert(xvalue_reduce)
    SX = np.fft.fft(xvalue_hilbert,n=NFFT)
    sx = SX[0:int(NFFT/2)]
    freq = fs/2*np.linspace(0,1,len(sx))
    # freq = freq[np.where(freq>0.1)]
    # sx = sx[np.where(freq>0.1)]
    return freq,sx 

def array_ifftx(t,sx):
    # fs = abs(f[1]-f[0])
    NFFT = 2*len(sx)
    # NFFT = int(2**(np.ceil(np.log2(len(f)))))
    # ts = 1/fs
    X = np.fft.ifft(sx,n=NFFT)
    # x = X[0:int(NFFT/2)]
    x = X[0:len(t)]
    # freq = freq[np.where(freq>0.1)]
    # sx = sx[np.where(freq>0.1)]
    return x

def ifftx(sxvalues,fvalues,pad):
    sxvalues = formatinput(sxvalues)
    fvalues = formatinput(fvalues)
    time = dict()
    x = dict()
    for name in list(sxvalues):
        dim1 = np.shape(sxvalues[name])
        dim2 = np.shape(fvalues[name])
        if len(dim1) == 1:
            fs = abs(fvalues[name][1]-fvalues[name][0])
            # NFFT = int(2**(np.ceil(np.log2(len(fvalues[name])))+pad))
            NFFT = 2*len(sxvalues[name])
            ts = 1/fs
            X = np.fft.ifft(sxvalues[name],n=NFFT)
            # x[name] = X[0:int(NFFT/2)]
            x[name] = X
            time[name] = ts/2*np.linspace(0,1,len(x[name]))
        if len(dim1) == 2:
            if len(dim2) == 1:
                fs = abs(fvalues[name][1]-fvalues[name][0])
                # NFFT = int(2**(np.ceil(np.log2(len(fvalues[name])))+pad))
                NFFT = 2*len(sxvalues[name])
                ts = 1/fs
                X = np.fft.ifft(sxvalues[name],n=NFFT,axis=0)
                # x[name] = X[0:int(NFFT/2)]
                x[name] = X
                time[name] = ts/2*np.linspace(0,1,len(x[name]))
            if len(dim2) == 2:
                fs = abs(fvalues[name][1,0]-fvalues[name][0,0])
                # NFFT = int(2**(np.ceil(np.log2(len(fvalues[name])))+pad))
                NFFT = 2*len(sxvalues[name])
                ts = 1/fs
                X = np.fft.ifft(sxvalues[name],n=NFFT,axis=0)
                # x[name] = X[0:int(NFFT/2)]
                x[name] = X
                time[name] = ts/2*np.linspace(0,1,len(x[name]))
    return x,time

def ifftx_omega(sxvalues,omega,pad=0):
    sxvalues = formatinput(sxvalues)
    omega = formatinput(omega)
    time = dict()
    x = dict()
    for name in list(sxvalues):
        dim1 = np.shape(sxvalues[name])
        dim2 = np.shape(omega[name])
        if len(dim1) == 1:
            fs = abs(omega[name][1]-omega[name][0])/2/np.pi
            NFFT = int(2**(np.ceil(np.log2(len(omega[name])))+pad))
            ts = 1/fs
            X = 2*np.pi*np.fft.ifft(sxvalues[name],n=NFFT)
            x[name] = X[0:int(NFFT/2)]
            time[name] = ts/2*np.linspace(0,1,len(x[name]))
        if len(dim1) == 2:
            if len(dim2) == 1:
                fs = abs(omega[name][1]-omega[name][0])/2/np.pi
                NFFT = int(2**(np.ceil(np.log2(len(omega[name])))+pad))
                ts = 1/fs
                X = 2*np.pi*np.fft.ifft(sxvalues[name],n=NFFT,axis=0)
                x[name] = X[0:int(NFFT/2)]
                time[name] = ts/2*np.linspace(0,1,len(x[name]))
            if len(dim2) == 2:
                fs = abs(omega[name][1,0]-omega[name][0,0])/2/np.pi
                NFFT = int(2**(np.ceil(np.log2(len(omega[name])))+pad))
                ts = 1/fs
                X = 2*np.pi*np.fft.ifft(sxvalues[name],n=NFFT,axis=0)
                x[name] = X[0:int(NFFT/2)]
                time[name] = ts/2*np.linspace(0,1,len(x[name]))
    return x,time

def dict_getabs(sxvalues):
    sx_out = {}
    names = list(sxvalues)
    for name in names:
        sx_out[name] = abs(sxvalues[name])
    return sx_out

def dict_getreal(sxvalues):
    names = list(sxvalues)
    for name in names:
        sxvalues[name] = np.real(sxvalues[name])
    return sxvalues

def array_FWHM(f, E):
    '''
    

    Parameters
    ----------
    freq : TYPE 1D numpy array
        DESCRIPTION. frequency
    E : TYPE
        DESCRIPTION. frequency domain field amplitude

    Returns
    -------
    FW : TYPE float
        DESCRIPTION. Full width half maximum

    '''
    hmax = max(E)/2
    idx0 = np.where(E >= hmax/2)[0][0]
    idx1 = np.where(E >= hmax/2)[0][-1]
    FW = abs(f[idx1]-f[idx0])
    return FW
        
def dict_key_get_T(key):
    key_strs = key.split('_')
    picked_str = ''
    for key_str in key_strs:
        if 'K' in key_str:
            picked_str = key_str
    if picked_str == '':
        num = 0
    else:
        num = ''
        for char in picked_str:
            if char.isdigit():
                num = num + char
        num = int(num)
    return num

def string_get_F(key):
    key_strs = key.split('_')
    picked_str = ''
    for key_str in key_strs:
        if 'F' in key_str:
            picked_str = key_str
            continue
    if picked_str == '':
        num = 0
    else:
        num = ''
        for char in picked_str:
            if char.isdigit():
                num = num + char
        num = int(num)
    return num

def dict_key_get_B(key):
    key_strs = key.split('_')
    picked_str = ''
    for key_str in key_strs:
        if 'T' in key_str:
            picked_str = key_str
    if picked_str == '':
        num = 0
    else:
        num = ''
        for char in picked_str:
            if char.isdigit():
                num = num + char
        num = int(num)
    return num
    

def dict_total_field(sx,sy):
    sall = dict() 
    for key in sx.keys():
        if key in list(sy):
            sall[key] = np.sqrt(sx[key]**2+sy[key]**2)
    return sall

def spec_sum(xvalues):
    xvalues = formatinput(xvalues)
    xsum = dict()
    for name in list(xvalues):
        xsum['SO'] = xsum['SO']+xvalues[name]
    return xsum

def array_getderivative(xvalues,yvalues):
    xdata = xvalues
    ydata = yvalues
    deri_store = np.zeros(len(xdata))
    for i in range(0,len(xdata)):
        if i == 0:
            deri_store[i] = (ydata[i+1]-ydata[i])/(xdata[i+1]-xdata[i])
        else:
            if i == len(xdata)-1:
                deri_store[i] = (ydata[i]-ydata[i-1])/(xdata[i]-xdata[i-1])
            else:
                deri_store[i] = (ydata[i+1]-ydata[i-1])/(xdata[i+1]-xdata[i-1])
    return deri_store
 
def dict_normalize(x):
    x = formatinput(x)
    xnew = dict()
    for name in x:
        xnew[name] = x[name]/max(x[name])
    return xnew
                         
def dict_square(x):
    x_sq = dict()
    for key,value in x.items():
        x_sq[key] = x[key]**2
    return x_sq



def find_list_common( str_list):
    str_list = list(str_list)
    com_str = ''
    if len(str_list) == 1:
        com_str = str_list[0]
    else:
        for i in range(0, len(str_list)-1):
            str1 = str_list[i]
            str2 = str_list[i+1]
            str1_s = str1.split('_')
            str2_s = str2.split('_')
            for chars in str1_s:
                if chars in str2_s:
                    if chars not in com_str:
                        com_str = com_str + chars + '_'
    return com_str

def find_str_common( str1, str2):
    com_str = ''
    str1_char = str1.split('_')
    str2_char = str2.split('_')
    for char in str1_char:
        if char in str2_char:
            com_str = com_str + char + '_'
    return com_str

def save_data(path, t, x, y = None):
    if y == None:
        if len(t) == len(x):
            y = np.zeros_like(t)
            spec = np.stack((t,x,y), axis = 1)
            np.savetxt(path, spec)
        else:
            raise Exception('dimension of time and amplitude does not match')
    else:
        if len(t) == len(x) and len(t) == len(y):
            spec = np.stack((t,x,y), axis = 1)
            np.savetxt(path, spec)
        else:
            raise Exception('Dimension of time and amplitude does not match')
         
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