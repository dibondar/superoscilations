# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 13:54:48 2023

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d


# t = np.linspace(-10, 10, 100) 
t = np.arange(-10,10, step = 0.05)
E = np.sin(5 * t)*np.exp(-(t)**2/5**2)

# Local freq estimation
def locfreq2(E):
    
    return np.sqrt(np.abs(
        np.gradient(np.gradient(E, t[1] - t[0]), t[1] - t[0]) / E
    ))
    #return np.sqrt(np.abs(
    #    gaussian_filter1d(E, 1, order=2)/ E    
    #))
    
locfreq = locfreq2(E)
plt.subplot(121)
plt.plot(t, E, label="FIeld")
plt.legend()

plt.subplot(122)
plt.plot(t[5:-5], locfreq2(E)[5:-5], )
plt.show()