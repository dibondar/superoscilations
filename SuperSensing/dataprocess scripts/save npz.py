# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 16:27:21 2023

@author: Admin
"""

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
        




path = 'C:/Users/Admin/Dropbox/Data/SPEC SO/2022/Dec/12-19-2022'
all_fnames = get_fnames(path)

t, x, y = load_single_spec(path, 'Glucose_SO')
np.savez(path+'/'+'Glucose_SO', time = t, amplitude = x)
t, x, y = load_single_spec(path, 'LGA_SO')
np.savez(path+'/'+'LGA_SO', time = t, amplitude = x)


for i in range(5,9):
    t, x, y = load_single_spec(path, 'Glucose_F0{}'.format(i))
    np.savez(path+'/'+'Glucose_F0{}'.format(i), time = t, amplitude = x)
    npz_save = np.load(path+'/'+'Glucose_F0{}.npz'.format(i))

    
for i in range(5,9):
    t, x, y = load_single_spec(path, 'LGA_F0{}'.format(i))
    np.savez(path+'/'+'LGA_F0{}'.format(i), time = t, amplitude = x)
    npz_save = np.load(path+'/'+'LGA_F0{}.npz'.format(i))
    
    
    
    
    