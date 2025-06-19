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


# def lactose_trans(freq, Ef, omega0, omega_p=0.069*2*pi*1e12, gamma=0.4*1e12, eps_inf=2.23):
#     c = 3e8
#     d = 1.10*1e-3
#     omega = freq*2*pi*1e12
#     omega0 = omega0*2*pi*1e12
#     omega0_2 = 5*2*pi*1e12
#     omega_p2 = 2*2*pi*1e12
#     gamma2 = 100*1e12
#     eps_lt = (eps_inf + omega_p ** 2 / (omega0 ** 2 - (omega) ** 2 - 1j * gamma * omega)
#               + omega_p2 ** 2 / (omega0_2 ** 2 - (omega) ** 2 - 1j * gamma2 * omega))
#     N = np.conj(np.sqrt(eps_lt))
#     E_out = Ef*2*N*2/(1+N)**2*np.exp(-1j*d*omega*(N-1)/c)*1/(1-(N-1)**2/(N+1)**2*np.exp(-2j*d*omega*N/c))
#     # a = -2/d*np.log(abs(E_out)/abs(Ef)*(1+N)**2/4/N)
#     k = np.imag(N)
#     n = np.real(N)
#     a = -4*pi*freq*k/c*1e12
#     # plt.figure()
#     # plt.plot(freq, Ef)
#     # plt.xlim((0,8))
#     # plt.figure()
#     plt.plot(freq, a)
#     plt.xlim((0,1.6))
#     return E_out

# def lga_trans(freq, Ef, omega0, omega_p=0.066*2*pi*1e12, gamma=0.4*1e12, eps_inf=2.04):
#     c = 3e8
#     d = 1.14*1e-3
#     omega = freq*2*pi*1e12
#     omega0 = omega0*2*pi*1e12
#     omega_p2 = 2.7*2*pi*1e12
#     omega0_2 = 5*2*pi*1e12
#     gamma2 = 40*1e12
#     eps_lt = (eps_inf + omega_p ** 2 / (omega0 ** 2 - (omega) ** 2 - 1j * gamma * omega)
#               + omega_p2 ** 2 / (omega0_2 ** 2 - (omega) ** 2 - 1j * gamma2 * omega))
#     N = np.conj(np.sqrt(eps_lt))
#     E_out = Ef*2*N*2/(1+N)**2*np.exp(-1j*d*omega*(N-1)/c)*1/(1-(N-1)**2/(N+1)**2*np.exp(-2j*d*omega*N/c))
#     # a = -2/d*np.log(abs(E_out)/abs(Ef)*(1+N)**2/4/N)
#     n = np.real(N)
#     k = np.imag(N)
#     a = -4*pi*freq*k/c*1e12
#     # plt.figure()
#     # plt.plot(freq, Ef)
#     # plt.xlim((0,8))
#     # plt.figure()
#     plt.plot(freq, a)
#     plt.xlim((0,1.6))
#     return E_out




bf = dp.Basic_functions()
# path = '/home/ppsapphire/Dropbox/SO project/artificial pulses'
path = 'C:/Users/Admin/Dropbox/SO project/compare low and high resonance for best J'
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


# Et.append(hilbert(np.sin(2*pi*0.5*time)*np.exp(-(time+10)**2/sigma**2)*1))
# Et.append(hilbert(np.sin(2*pi*0.6*time)*np.exp(-(time+10)**2/sigma**2)*1))
# Et.append(hilbert(np.sin(2*pi*0.7*time)*np.exp(-(time+10)**2/sigma**2)*1))
# Et.append(hilbert(np.sin(2*pi*0.8*time)*np.exp(-(time+10)**2/sigma**2)*1))
freqs = np.arange(0.5, 0.81, step = 0.1)

for f in freqs:
    Et.append(np.cos(2*pi*f*time)*np.exp(-(time)**2/sigma**2)*1)
# Et.append(np.cos(2*pi*0.5*time)*np.exp(-(time)**2/sigma**2)*1)
# Et.append(np.sin(2*pi*0.53*time)*np.exp(-(time)**2/sigma**2)*1)
# Et.append(np.sin(2*pi*0.53*time)*np.exp(-(time)**2/sigma**2)*1)
# Et.append(np.cos(2*pi*0.6*time)*np.exp(-(time)**2/sigma**2)*1)
# Et.append(np.sin(2*pi*0.65*time)*np.exp(-(time)**2/sigma**2)*1)
# Et.append(np.cos(2*pi*0.7*time)*np.exp(-(time)**2/sigma**2)*1)
# Et.append(np.sin(2*pi*0.75*time)*np.exp(-(time)**2/sigma**2)*1)
# Et.append(np.cos(2*pi*0.8*time)*np.exp(-(time)**2/sigma**2)*1)




# Et.append(np.sin(2*pi*0.5*time))
# Et.append(-np.sin(2*pi*0.6*time))
# Et.append(np.sin(2*pi*0.7*time))
# Et.append(-np.sin(2*pi*0.8*time))


plt.figure()
plt.plot(time,Et[0])

Et_tx = []
Et_ty = []
Ets_tx = []
Ets_ty = []
Et2_tx = []
Et2_ty = []
Et2s_tx = []
Et2s_ty = []
    
ampl_max = []

#calculate transmission coefficient of InSb for x and y component
# freq, Ef = bf.array_fftx(Et[0], time)
# Tx, Ty = np.nan_to_num(FaradayRotate(freq,0.125))
# Tx2, Ty2 = np.nan_to_num(FaradayRotate(freq,0.128))

# calculate wave transmitted through InSb
# for i in range(0,4):
#     freq, Ef = bf.array_fftx(Et[i], time)
#     Et_tx.append(np.real(bf.array_ifftx(time, Ef*Tx)[idx1:idx2]))
#     Et_long_tx.append(np.real(bf.array_ifftx(time, Ef*Tx)[0:len(time)]))
#     Et_ty.append(np.real(bf.array_ifftx(time, Ef*Ty)[idx1:idx2]))
#     Et_long_ty.append(np.real(bf.array_ifftx(time, Ef*Ty)[0:len(time)]))
#     Et2_tx.append(np.real(bf.array_ifftx(time, Ef*Tx2)[0:len(time)]))
#     Et2_ty.append(np.real(bf.array_ifftx(time, Ef*Ty2)[0:len(time)]))
#     ampl_max.append(max(np.real(bf.array_ifftx(time, Ef*Ty)[idx1:idx2])))

for i in range(0,len(freqs)):
    freq, Ef = bf.array_fftx_hilbert(Et[i], time)
    # Et_ty.append(np.real(bf.array_ifftx(time, lactose_trans(freq, Ef, omega0=inf))[idx1:idx2]))
    # Et_tx.append(np.zeros(len(Et_ty[i])))
    Et_ty.append(np.real(bf.array_ifftx(time, lga_trans(freq, Ef, omega0=1.2))))
    Et_tx.append(np.zeros(len(Et_ty[i])))
    Ets_ty.append(np.real(bf.array_ifftx(time, lga_trans(freq, Ef, omega0=1.2)))[idx1:idx2])
    Ets_tx.append(np.zeros(len(Ets_ty[i])))
    # Et2_ty.append(np.real(bf.array_ifftx(time, lactose_trans(freq, Ef, omega0=1.21*1e12))[0:len(time)]))
    Et2_ty.append(np.real(bf.array_ifftx(time, lga_trans(freq, Ef, omega0=1.4))))
    Et2_tx.append(np.zeros(len(Et2_ty[i])))
    Et2s_ty.append(np.real(bf.array_ifftx(time, lga_trans(freq, Ef, omega0=1.4)))[idx1:idx2])
    Et2s_tx.append(np.zeros(len(Et2s_ty[i])))
    
    # Et_long_ty.append(np.real(bf.array_ifftx(time, lactose_trans(freq, Ef, omega0=inf))))
    ampl_max.append(max(Et_ty[i]))
    # Et_long_tx.append(np.zeros(len(Et_long_ty[i])))
plt.figure()
plt.plot(time, Et[0], time, Et_ty[0], time, Et2_ty[0])
    # plt.xlabel('time (ps)')
    # plt.ylabel('amplitude (a.u.)')
    
ampl_max = np.array(ampl_max)
ampl_ult = max(ampl_max)

# adjust amplitude of the individual pulses to give better superoscillation
# ampl_factor = np.array([1,1,1,1])
# for i in range(0,4):
#     coef = ampl_ult/ampl_max[i]
#     Et_tx[i] = Et_tx[i]*coef*ampl_factor[i]
#     Et_ty[i] = Et_ty[i]*coef*ampl_factor[i]
#     Ets_tx[i] = Et1_tx[i]*coef*ampl_factor[i]
#     Ets_ty[i] = Et1_ty[i]*coef*ampl_factor[i]
#     Et2_tx[i] = Et2_tx[i]*coef*ampl_factor[i]
#     Et2_ty[i] = Et2_ty[i]*coef*ampl_factor[i]
    
# plt.figure()
# plt.plot(time_s,Et_ty[0],time_s,Et_ty[1],time_s,Et_ty[2],time_s,Et_ty[3])
# plt.xlabel('time (ps)')
# plt.ylabel('amplitude (a.u.)')

# plt.figure()
# plt.plot(time_s,Et_tx[0],time_s,Et_tx[1],time_s,Et_tx[2],time_s,Et_tx[3])
# plt.xlabel('time (ps)')
# plt.ylabel('amplitude (a.u.)')

# plt.figure()
# plt.plot(time,Et_long_ty[0],time,Et_long_ty[1],time,Et_long_ty[2],time,Et_long_ty[3])
# plt.xlabel('time (ps)')
# plt.ylabel('amplitude (a.u.)')

for i in range(0,len(freqs)):
    spec0 = np.stack((time, Et[i], Et_tx[i]), axis = 1)
    spec = np.stack((time, Et_ty[i], Et_tx[i]), axis=1)
    spec2 = np.stack((time, Et2_ty[i], Et2_tx[i]), axis=1)
    spec3 = np.stack((time_s, Ets_ty[i], Ets_tx[i]), axis=1)
    spec4 = np.stack((time_s, Et2s_ty[i], Et2s_tx[i]), axis = 1)
    # spec = spec.transpose(1,0)
    flabel = str(int(freqs[i]*100))
    for j in range(0,3):
        
        # fname0 = path+'/incident_F0'+str(i+5)+'_'+str(j+1)+'.dat'
        filename = path+'/H_N0_F'+str(flabel)+'_'+str(j+1)+'.dat'
        filename2 = path+'/H_N2_F'+str(flabel)+'_'+str(j+1)+'.dat'
        filename3 = path+'/H_N0S_F'+str(flabel)+'_'+str(j+1)+'.dat'
        filename4 = path+'/H_N2S_F'+str(flabel)+'_'+str(j+1)+'.dat'
        # np.savetxt(fname0, spec0)
        np.savetxt(filename, spec)
        np.savetxt(filename2, spec2)
        np.savetxt(filename3, spec3)
        np.savetxt(filename4, spec4)


