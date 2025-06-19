# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 09:58:25 2023

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

def FaradayRotate(freq,B):
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
    T_x = (T_cra+T_cri)/np.sqrt(2)
    T_y = (T_cra-T_cri)/np.sqrt(2)/1j
        
    return T_x, T_y


def lactose_trans(freq, Ef, omega0, omega_p=0.05*2*pi*1e12, gamma=0.1*1e12, eps_inf=2.3):
    c = 3e8
    d = 1.14*1e-3
    omega = freq*2*pi*1e12
    omega0 = omega0*2*pi*1e12
    eps_lt = eps_inf + omega_p ** 2 / (omega0 ** 2 - (omega) ** 2 - 1j * gamma * omega)
    N = np.conj(np.sqrt(eps_lt))
    E_out = Ef*2*N*2/(1+N)**2*np.exp(-1j*d*omega*(N-1)/c)*1/(1-(N-1)**2/(N+1)**2*np.exp(-2j*d*omega*N/c))
    # a = -2/d*np.log(abs(E_out)/abs(Ef)*(1+N)**2/4/N)
    k = np.imag(N)
    a = -4*pi*freq*k/c*1e12
    # plt.figure()
    # plt.plot(freq, Ef)
    # plt.xlim((0,8))
    # plt.figure()
    plt.plot(freq, a)
    plt.xlim((0,1.6))
    return E_out

def lga_trans(freq, Ef, omega0, omega_p=0.05*2*pi*1e12, gamma=0.1*1e12, eps_inf=2.3):
    c = 3e8
    d = 1.14*1e-3
    omega = freq*2*pi*1e12
    omega0 = omega0*2*pi*1e12
    eps_lt = eps_inf + omega_p ** 2 / (omega0 ** 2 - (omega) ** 2 - 1j * gamma * omega)
    N = np.conj(np.sqrt(eps_lt))
    E_out = Ef*2*N*2/(1+N)**2*np.exp(-1j*d*omega*(N-1)/c)*1/(1-(N-1)**2/(N+1)**2*np.exp(-2j*d*omega*N/c))
    # a = -2/d*np.log(abs(E_out)/abs(Ef)*(1+N)**2/4/N)
    k = np.imag(N)
    a = -4*pi*freq*k/c*1e12
    # plt.figure()
    # plt.plot(freq, Ef)
    # plt.xlim((0,8))
    # plt.figure()
    plt.plot(freq, a)
    plt.xlim((0,1.6))
    return E_out





bf = dp.Basic_functions()
# path = '/home/ppsapphire/Dropbox/SO project/artificial pulses'
path = 'C:/Users/Admin/Dropbox/SO project/artificial pulses for new J'
# path = '/Users/ppsap/Dropbox/SO project/artificial pulses'
path_spec = 'C:/Users/Admin/Dropbox/SO project/artificial pulses'


'''read in pulses'''

filename = 'export_Glucose_5per_295K_ref_1.dat'
spec1 = np.loadtxt(path_spec + '/' + filename)
# spec2 = np.loadtxt(path + '/' + filename2)

t1 = spec1[:,0]
t_shift = t1 - 551
# t2 = spec2[:,0]
x1 = spec1[:,1]
# x2 = spec2[:,1]
plt.figure()
plt.plot(t1, x1)

f1, sx1 = bf.array_fftx_hilbert(x1, t1, pad = 2)
# f2, sx2 = bf.array_fftx_hilbert(x2, t2, pad = 2)

plt.figure()
x1_T = np.real(bf.array_ifftx(t1, lga_trans(f1, sx1, omega0=1.2)))
x2_T = np.real(bf.array_ifftx(t1, lactose_trans(f1, sx1, omega0=1.43)))

zero_filler = np.zeros_like(t1)

spec1_new = np.stack((t_shift, x1_T, zero_filler), axis = 1)
spec2_new = np.stack((t_shift, x2_T, zero_filler), axis = 1)

np.savetxt(path + '/' + 'broad_LGA_1.dat', spec1_new)
np.savetxt(path + '/' + 'broad_Glucose_1.dat', spec2_new)

plt.figure()
plt.plot(t1, x1, t1, x1_T, t1, x2_T)






