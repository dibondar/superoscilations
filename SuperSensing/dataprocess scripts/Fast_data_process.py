# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 14:04:00 2023

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
# from numba import njit, jit
from numpy.fft import fft


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
    dim_y = 1
    for i in range(0, len(flist)):
        fnum = str(i+1)
        try:
            data = np.loadtxt(folder+'/'+fname+'_'+fnum+ftype)
        except OSError:
            continue
        else:
            data = np.loadtxt(folder+'/'+fname+'_'+fnum+ftype)
            if i == 0:
                t_store = data[:, 0]
                x_store = data[:, 1]
                try:
                    y_store = data[:, 2]
                except:
                    dim_y = 0
                else:
                    y_store = data[:, 2]
            else:
                if dim_y == 0:
                    t_new = data[:, 0]
                    x_new = data[:, 1]
                    t_store = np.colume_stack((t_new, t_store))
                    x_store = np.column_stack((x_new, x_store))
                else:
                    t_new = data[:, 0]
                    x_new = data[:, 1]
                    y_new = data[:, 2]
                    t_store = np.column_stack((t_new, t_store))
                    x_store = np.column_stack((x_new, x_store))
                    y_store = np.column_stack((y_new, y_store))
    return t_store, x_store, y_store

def str_get_D(string):
    str_seg_all = string.split('_')
    pick_str = ''
    for str_seg in str_seg_all:
        if 'D' in str_seg:
            pick_str = str_seg
            continue
    num = ''
    for char in pick_str:
        if char.isdigit():
            num = num + char
    try:
        num = int(num)
    except ValueError:
        num = 0
    else:
        num = int(num)
    return num




def str_get_Delay(string):
    str_seg_all = string.split('_')
    pick_str = ''
    for str_seg in str_seg_all:
        if 'Delay' in str_seg:
            pick_str = str_seg
            continue
    num = ''
    for char in pick_str:
        if char.isdigit():
            num = num + char
    try:
        num = int(num)
    except ValueError:
        num = 0
    else:
        num = int(num)
    return num

# @njit
def interferometric_vis(sx):
    # sx_normal = abs(sx)/max(abs(sx))
    sx_normal = np.absolute(sx)
    i_max = np.max(sx_normal**2)
    i_min = np.min(sx_normal**2)
    ivis = (i_max-i_min)/(i_max+i_min)
    return ivis


def average_spec_noref(t, x):
    x_ave = {}
    t_ave = {}
    for key in list(x):
        try:
            t_ave[key] = t[key][:,0]
        except IndexError:
            t_ave[key] = t[key]
            x_ave[key] = x[key]
        else:
            t_ave[key] = t[key][:,0]
            x_ave[key] = np.average(x[key], axis = 1)
    return t_ave, x_ave

# @njit
def generate_gauss_cos(t, f, sigma = 1):
    return np.cos(2*pi*f*t)*np.exp(-(t)**2/sigma**2)


def average_spec_samref(x, t, nsam, nref):
    samx = dict()
    samt = dict()
    refx = dict()
    reft = dict()
    comp1 = dict()
    comp2 = dict()
    for key in list(x):
        key = key
        x_temp = x[key]
        t_temp = t[key]
        t_temp = t_temp
        samx_sum = 0
        refx_sum = 0
        n_cycle = round(np.size(x_temp, axis = 1)/(nsam+nref))
        for i in range(np.size(x_temp, axis=1)):
            if i % (nsam+nref) < nsam:
                samx_sum = x_temp[:,i] + samx_sum
            else:
                refx_sum = x_temp[:,i] + refx_sum
        samx[key] = samx_sum/n_cycle
        samt[key] = t_temp[:,0]
        reft[key] = t_temp[:,nref]
        refx[key] = refx_sum/n_cycle
        comp1[key] = samx[key]+refx[key]
        comp2[key] = samx[key]-refx[key]
    return samx,samt,refx,reft,comp1,comp2
        
    


def fftx(t, x, pad=2):
    ts = abs(t[1]-t[0])
    NFFT = int(2**(np.ceil(np.log2(len(t)))+pad))
    fs = 1/ts
    x_reduce = x-np.average(x)
    SX = np.fft.fft(x_reduce,n=NFFT)
    sx = SX[0:int(NFFT/2)]
    freq = fs/2*np.linspace(0,1,len(sx))
    return freq, sx 


def ifftx(t,sx):
    NFFT = 2*len(sx)
    X = np.fft.ifft(sx,n=NFFT)
    x = X[0:len(t)]
    return x


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


def fftx_stft(t, x, pad = 2):
    ts = abs(t[1]-t[0])
    NFFT = int(2**(np.ceil(np.log2(len(t)))+pad))
    fs = 1/ts
    x_normal = x-np.average(x)
    f_stft, t_stft, zx = stft(x_normal, fs, nperseg = len(x)/2, noverlap=len(x)/2-1,
                                        window = ('gaussian',14), nfft = NFFT)
    return f_stft, t_stft, zx

def localf_stft(t, x, pad = 2):
    ts = abs(t[1]-t[0])
    NFFT = int(2**(np.ceil(np.log2(len(t)))+pad))
    fs = 1/ts
    x_normal = x-np.average(x)
    f_stft, t_stft, zxx = stft(x_normal, fs, nperseg = len(x)/2, noverlap=len(x)/2-1,
                                        window = ('gaussian',14), nfft = NFFT)
    nrow, ncol = np.shape(zxx)
    ampl = abs(zxx)
    f_ave = np.zeros(ncol)
    for i in range(0, ncol):
        ampl_sum = sum(ampl[:,i])
        f_ave[i] = (sum(f_stft*ampl[:,i])/ampl_sum) 
    return t_stft, f_ave


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


def spec_shift(t, x, t_shift):
    step = abs(t[1]-t[0])
    move = int(np.floor(t_shift/step))
    if t_shift > 0:
        x_new = np.concatenate((np.zeros(move), x[0:-move]), axis = 0)
    if t_shift < 0:
        x_new = np.concatenate((x[move:], np.zeros(move)))
    return x_new
    
    


if __name__ == '__main__':
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