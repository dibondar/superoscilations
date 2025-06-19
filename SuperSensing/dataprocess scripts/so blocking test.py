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
from numba import njit
import plottool_v6 as pt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
import re
# pad = np.zeros(100)
# t = np.arange(0, 1.2, step = 0.05)
# x = np.sin(2*pi*2*t)
# f, sx = dp.Basic_functions().array_fftx(x, t)
# fig, ax = plt.subplots(1,1,figsize = (9,6), dpi=300)
# ax.plot(f, abs(sx))

def get_data_fnames(path, filetype='.dat'):
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
path = 'C:/Users/Admin/Dropbox/Data/SPEC SO/2023/Jun/06-22-2023'
# fnames = get_data_fnames(path)
fnames = ['Block_SO']
# t_all, x_all, y_all = dp.Basic_functions.load_single_spec(path, 'Block_SO')
distance = np.linspace(6, 11, num = 51)
integral_ci_all = {}
integral_all = {}
gradient_all = {}
distance_all = {}
for fname in fnames:
    t_all, x_all, y_all = load_single_spec(path, fname)
    ndatapoints, nscan = np.shape(t_all)
    integral_so = []
    integral_ci = []
    t0 = t_all[:,0]
    x0 = x_all[:,0]
    # x_abs = abs(x0)
    t_m = t0[np.where(x0 == max(x0))[0][0]]
    for i in range(0, nscan):
        t = t_all[:, i]
        x = x_all[:, i]
        x_intense = x**2
        x_so = x_intense[(t >= -0.6) & (t <= 0.6)]
        # idx_so0 = np.where(t>=-0.6)[0][0]
        # idx_so1 = np.where(t<=0.6)[0][-1]
        x_ci = x_intense[(t >= t_m-0.6) & (t <= t_m+0.6)]
        integral_so.append(simps(x_so, dx = t[1]-t[0]))
        integral_ci.append(simps(x_ci, dx = t[1]-t[0]))
    integral_ci = np.array(integral_ci)
    integral_so = np.array(integral_so)
    ci_gradient = np.gradient(integral_ci, distance)
    so_gradient = np.gradient(integral_so, distance)
    ci_gradient = ci_gradient/max(abs(ci_gradient))
    so_gradient = so_gradient/max(abs(so_gradient))
    integral_all[fname+'_ci'] = integral_ci/max(integral_ci)
    integral_all[fname+'_so'] = integral_so/max(integral_so)
    integral_all[fname+'_ci_gradient'] = ci_gradient
    integral_all[fname+'_so_gradient'] = so_gradient
    distance_all[fname+'_ci'] = distance
    distance_all[fname+'_so'] = distance

root = tk.Tk()
window = tk.Frame(root)
window.grid()
plotgo = pt.Plottool(window, distance_all, integral_all)     
root.mainloop()  
# integral_so = np.array(integral_so)

# integral_so = integral_so/max(integral_so)


# plt.plot(distance, integral_so, label = 'Superoscillation')
# plt.plot(distance, integral_ci, label = 'constructive interference')
# plt.legend()
# plt.xlabel('Bloack travel (mm)')
# plt.ylabel('$\int(E^2)dt$')
# arr = np.array(([1,2,3],[1,2,3]))
# ncol = arr.size


# t ,x,y = dp.Basic_functions.load_single_spec(path, 'PartialBlock_2nd_SO')
# flist = glob.glob(path+'/'+'PartialBlock_2nd_SO'+'_*[0-9]'+'.dat')
# fnum = str(1)
# data = np.loadtxt(path+'/'+'PartialBlock_2nd_SO'+'_'+fnum+'.dat')


