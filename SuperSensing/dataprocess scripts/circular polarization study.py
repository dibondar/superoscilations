# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 16:28:03 2023

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import cos
from numpy import sin
from numpy import pi
import dataprocessor as dp
import Fast_data_process as fdp

def Faraday_transmission(freq, B, t, Ex, Ey):
    '''
    

    Parameters
    ----------
    freq : 1d array
        DESCRIPTION. frequency domain obtained from Fourier transform input waves
    B : TYPE float64
        DESCRIPTION. External magnetic field strength in T

    Returns
    -------
    T_x : 1d array
        DESCRIPTION. complex transmission coefficient in x direction
    T_y : 1d array
        DESCRIPTION. complex transmission coefficient in y direction

    '''
    
    d=5*1e-4 #sample thickness
    w=2*pi*freq*1e12 #angular frequency 
    #set parameters for calculation
    n_air=1 # index of air
    e=1.6e-19 #electron charge
    c=3e8 # speed of light
    v_e=0.86*1e12 # electron scattering rate
    m_e=0.018*9.1e-31 # mass of electron
    w_t=2*pi*5.90e12 # transverse phonon frequency
    w_l=2*pi*5.54e12 # longitudinal phonon frequency
    epi_b=13.82 # background permittivity
    w_ce=e*B/(m_e) # cyclotron frequency
    w_e=0.28*2*pi*1e12*np.sqrt(epi_b) #plasma frequency is set to constant 0.28
    # w_h=(4*pi*Nh*e^2/m_h*1);
    v_ph=2*pi*1e12 # phonon damping rate
    
    fx, Ex_f = fdp.fftx(t, Ex)
    fy, Ey_f = fdp.fftx(t, Ey)
    Ef_cra = Ex_f+1j*Ey_f
    Ef_cri = Ex_f-1j*Ey_f
   
    epi_ph = epi_b*((w_t**2-w_l**2)/(w_t**2-w**2-1j*v_ph*w))

    ncra=np.sqrt(epi_b-w_e**2/(w*(w+1j*v_e-w_ce))+epi_ph)
    ncri=np.sqrt(epi_b-w_e**2/(w*(w+1j*v_e+w_ce))+epi_ph)

    t_cra_as=2*n_air/(n_air+ncra) #calculate Fresnel transmission coefficient
    t_cra_sa=2*ncra/(n_air+ncra)
    t_cri_as=2*n_air/(n_air+ncri)
    t_cri_sa=2*ncri/(n_air+ncri)
    # calculate complex transmission coefficents
    T_cra =t_cra_as*t_cra_sa*np.exp(1j*(ncra-1)*w*d/c)
    T_cri =(t_cri_as*t_cri_sa*np.exp(1j*(ncri-1)*w*d/c))
    Ef_cra_T = T_cra*Ef_cra
    Ef_cri_T = T_cri*Ef_cri
    Ef_x = (Ef_cri_T+Ef_cra_T)/np.sqrt(2)
    Ef_y = (Ef_cra_T-Ef_cri_T)/np.sqrt(2)/1j
    Et_x = fdp.ifftx(t, Ef_x)
    Et_y = fdp.ifftx(t, Ef_y)
    return Et_x, Et_y

path = 'C:/Users/Admin/Dropbox/SO project/simulation for circular polarized SO'

t = np.arange(-20,20,step = 0.05)
ts = np.arange(-5,5,step = 0.05)
# omega_all = np.array([0.5, 0.6, 0.7, 0.8])*1e12*2*pi
freq_all = np.array([0.5, 0.6, 0.7, 0.8])
# omega = 0.8*1e12*2*pi
# omega_1 = 0.1*1e12*2*pi
Ey_L = []
Ex = {}
Ey = {}
sigma = 20
Et_x = 0
Et_y = 0
# generate original pulses

for idx, freq in enumerate(freq_all, start = 5): 
    omega = freq*2*pi
    Ex_temp = cos(omega*t)*np.exp(-(t+0)**2/sigma**2)
    Ey_temp = np.zeros_like(Ex_temp)
    # Ex_temp = (cos(omega*t)-sin(omega*t))*np.exp(-(t+0)**2/sigma**2)
    # Ey_temp = (cos(omega*t)+sin(omega*t))*np.exp(-(t+0)**2/sigma**2)
    key = 'F0'+str(idx)
    # Ex[key] = Ex_temp
    # Ey[key] = Ey_temp
    Ex_T, Ey_T = Faraday_transmission(freq, 0.125, t, Ex_temp, Ey_temp)
    Ex_T2, Ey_T2 = Faraday_transmission(freq, 0.130, t, Ex_temp, Ey_temp)
    
    spec1 = np.stack((t, Ex_T, Ey_T), axis = 1)
    spec2 = np.stack((t, Ey_T, Ex_T), axis = 1)
    spec3 = np.stack((t, Ex_T2, Ey_T2), axis = 1)
    spec4 = np.stack((t, Ey_T2, Ex_T2), axis = 1)
    plt.figure()
    plt.plot(t, Ex_T, t, Ex_T2)
    plt.figure()
    plt.plot(t, Ey_T, t, Ey_T2)
    
    for i in range(1,4):
        np.savetxt(path+'/CX0_'+key+'_{}.dat'.format(i), spec1)
        np.savetxt(path+'/CY0_'+key+'_{}.dat'.format(i), spec2)
        np.savetxt(path+'/CX2_'+key+'_{}.dat'.format(i), spec3)
        np.savetxt(path+'/CY2_'+key+'_{}.dat'.format(i), spec4)
    
    # Ey_L.append(np.sqrt(2)*cos(omega*t))
    # plt.figure()
    # ax = plt.axes(projection = '3d')
    # ax.plot(Ex_temp, Ey_temp, t)

# for key in list(Ex):
#     Et_x = Et_x + Ex[key]
#     Et_y = Et_y + Ey[key]
    


# plt.figure()
# plt.plot(t, Et_x)
# plt.figure()
# plt.plot(t, Et_y)
    

  

