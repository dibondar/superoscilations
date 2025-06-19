# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 16:14:40 2023

@author: Admin
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 10:55:04 2022

@author: ppsapphire
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy.signal import hilbert
import dataprocessor as dp


bf = dp.Basic_functions()
# path = '/home/ppsapphire/Dropbox/SO project/artificial pulses'
path = 'C:/Users/Admin/Dropbox/SO project/artificial pulses for more FHs'
# path = '/Users/ppsap/Dropbox/SO project/artificial pulses for best ND'



Et = []
time = np.arange(-55, 55, step = 0.05)

#set small time window for SO phase calculation
idx1 = np.where(time>=-20)[0][0]
idx2 = np.where(time<=20)[0][-1]
time_s = time[idx1:idx2]


#generate decaying sin waves
sigma = 20
inf = float('inf')


for i in range(0,4):
    Et.append(np.cos(2*pi*0.6*time)*np.exp(-(time)**2/sigma**2)*1)


plt.figure()
plt.plot(time,Et[0])

Et_tx = []
Et_ty = []
Ets_tx = []
Ets_ty = []





for i in range(0,4):
    Et_ty.append(Et[i])
    Et_tx.append(np.zeros_like(Et[i]))
    Ets_ty.append(Et[i][idx1:idx2])
    Ets_tx.append(np.zeros_like(Et[i][idx1:idx2]))




for i in range(0,4):
    spec0 = np.stack((time, Et_ty[i], Et_tx[i]), axis = 1)
    spec1 = np.stack((time_s, Ets_ty[i], Ets_tx[i]), axis=1)

    # spec = spec.transpose(1,0)
    for j in range(0,3):
        # fname0 = path+'/incident_F0'+str(i+5)+'_'+str(j+1)+'.dat'
        filename0 = path+'/P_N'+str(i)+'_F06'+'_'+str(j+1)+'.dat'
        filename1 = path+'/P_SN'+str(i)+'_F06'+'_'+str(j+1)+'.dat'
        np.savetxt(filename0, spec0)
        np.savetxt(filename1, spec1)
     


