import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import tkinter as tk
import dataprocessor as dp
import Fast_data_process as fdp
import glob

import copy

from scipy import stats
from scipy.signal import find_peaks, hilbert, windows, stft
from scipy.integrate import simps
from scipy.interpolate import UnivariateSpline, InterpolatedUnivariateSpline
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize


from collections import Counter
from itertools import product
from dataclasses import dataclass

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk

from multiprocessing import Pool, cpu_count, Process
from tqdm.notebook import tqdm
import configparser as cfp





@dataclass
class Pulse:
    time: np.ndarray
    field: np.ndarray
    _individual_fields: np.ndarray = None
    peaks_time: np.ndarray = None
    time_range: np.ndarray = None
    half_period: float = 0
    interp_field: InterpolatedUnivariateSpline = None
    err_lower_bound: InterpolatedUnivariateSpline = None
    err_upper_bound: InterpolatedUnivariateSpline = None
 

class Calculate_SO_phase:
    def __init__(self, path, filenames_wf_base, filetype = '.dat', progress = None):
        self.pulses = {}
        self.max_ampl = []
        fnames = sorted(filenames_wf_base, key = dp.Basic_functions.string_get_F)
        #load data
        for name in fnames:
            # load all pulses for the given frequency
            self.filenames = glob.glob(path + name + '_*[0-9].dat')
            self.fnames_valid = []
            for j in range(0,len(self.filenames)):
                filenum = str(j+1)
                try:
                    data = np.loadtxt(path+'/'+name+'_'+filenum+filetype)
                except OSError:
                    continue
                else:
                    self.fnames_valid.append(path+'/'+name+'_'+filenum+filetype)
            self.all_data = [np.loadtxt(_).T[:-1] for _ in self.fnames_valid]
            self.times, self._individual_fields = zip(*self.all_data)
            self._individual_fields = np.array(self._individual_fields)
            time_mid = round((self.times[0][0]+self.times[0][-1])/2)
            time_diff = time_mid - 0
            self.time = self.times[0] - time_diff
            self.label = name
            self.field = np.mean(self._individual_fields, axis=0)
            # Substract the mean from experimental fields to compensate for parasitic offset
            self.field -= self.field.mean()
            # calculate the confidence interval
            confidence_level = 0.95
            self.sigma = stats.sem(self._individual_fields, axis=0)
            self.sigma = np.clip(self.sigma, 1e-20, self.sigma.max())
            self.err_lower_bound, self.err_upper_bound = stats.t.interval(confidence_level, 
                                                                          self._individual_fields.shape[0], 
                                                                          self.field, self.sigma,)
            abs_field = np.abs(self.field)
            self.max_ampl.append(abs_field.max())
            # Extract information for the combinatorial method     
            ampl_threshold = 0.8 * abs_field.max()
    
            indx = find_peaks(abs_field, height=ampl_threshold)[0]
            peaks_time = self.time[indx[1:-1]]
    
            self.half_period = Counter(np.diff(peaks_time)).most_common(1)[0][0]
        
            
                # Saving the data 
            self.pulses[self.label] = Pulse(
                time = self.time, 
                field = self.field,
                _individual_fields = self._individual_fields,
                peaks_time = peaks_time,
                time_range = self.time[indx[1]:indx[-2]],
                half_period = self.half_period,
                # first save confidence interval as arrays
                err_lower_bound = self.err_lower_bound,
                err_upper_bound = self.err_upper_bound,
            )
            
            self.fastest_field = self.field
        
        
        self.largest_freq = self.label

        #  Checking whether the time axis coincide
        #assert all(np.allclose(time, data.time) for data in pulses.values()), \
            #"This workbook cannot be used since the time data is not syncronized"

        # saving time step
        self.dtime = self.time[1] - self.time[0]

        self.max_ampl = max(self.max_ampl)
            
        

        # Normalazing fields and interpolating
        for data in self.pulses.values():
            data.field /= self.max_ampl
            
            data.interp_field = UnivariateSpline(
                data.time, 
                data.field,
                ext='zeros', 
                k=3, 
                s=0
            )
            
        
            # interpolate confidence interval
            data.err_lower_bound = UnivariateSpline(
                data.time, 
                data.err_lower_bound /self.max_ampl, 
                ext='zeros', 
                k=3, 
                s=0
            )
            
            data.err_upper_bound = UnivariateSpline(
                data.time, 
                data.err_upper_bound / self.max_ampl, 
                ext='zeros', 
                k=3, 
                s=0
            )
         
            
            #print(data.field)
            #print(list(pulses))
        self.observational_window = -0.6, 0.6
        self.half_period = self.pulses[self.largest_freq].half_period
        self.time_window, self.dx = np.linspace(-0.6, 0.6, 100, retstep=True)
        
        self.time_window_raw = self.time[(self.time >= -0.6) & (self.time <= 0.6)] 
        
        # get so quantify number for fastest oscillation
        self.fastest_field_short = self.fastest_field[(self.time >= -0.6) & (self.time <= 0.6)]
        self.fast_freq, self.fastest_field_fd = dp.Basic_functions().array_fftx(self.fastest_field_short, self.time_window_raw)
        self.fast_freq_max = self.fast_freq[np.where(abs(self.fastest_field_fd) == max(abs(self.fastest_field_fd)))][0]
        self.fastest_FWHM = dp.Basic_functions.array_FWHM(self.fast_freq, abs(self.fastest_field_fd))
        self.so_qnum = self.fast_freq_max/self.fastest_FWHM
        #generate initial guesses
        np.random.seed(3112022)

        self.bounds = [(_.peaks_time.min(), _.peaks_time.max()) for _ in self.pulses.values()]
        
        #bounds_mid = [(_.peaks_time.min()+ _.peaks_time.max())/2 for _ in pulses.values()]
        #print(bounds_mid)
        self.rand_initial_guess = np.array([
            np.random.uniform(_[0], _[1], 400) for _ in self.bounds
        ])
        # self.rand_initial_guess_test = self.rand_initial_guess
        # for i in range(0,len(self.rand_initial_guess)):
        #     self.rand_initial_guess[i,0] = 0
        
        #print(rand_initial_guess)
        self.rand_initial_guess = set(tuple(_) for _ in self.rand_initial_guess.T)
        
        pool = Pool(processes = int(cpu_count()))
        
        self.gradient_descent_results = set(pool.map(self.local_minimization, self.rand_initial_guess))
            
        # self.gradient_descent_results = set(self.local_minimization(_) for _ in (self.rand_initial_guess))
        
       
        self.gradient_descent_results = sorted(self.gradient_descent_results)
        self.intensity_without_ampl_modulation, self.all_time_delays = zip(*self.gradient_descent_results)
        
        self.all_time_delays = set(tuple(_) for _ in self.all_time_delays)
        self.so_quant_results = []
        for time_delays in self.all_time_delays:
            self.so_quant_results.append(self.quant_so(time_delays))
        self.so_quant_results = sorted(self.so_quant_results, reverse = True)
        self.so_qnums, self.all_time_delays = zip(*self.so_quant_results)
        

    def get_combined_field(self, time_delays, time_window):
        return sum(_.interp_field(time_window - delay) for delay, _ in zip(time_delays, self.pulses.values()))

    def get_combined_field_test(self, time_delays, time_window):
        new_delay = (0, time_delays[1], time_delays[2], time_delays[3])
        return sum(_.interp_field(time_window - delay) for delay, _ in zip(new_delay, self.pulses.values()))
        
    
            

    def inegral_without_ampl_modulation(self,time_delays):
        E_c = self.get_combined_field(time_delays, self.time_window)
        freq, E_f = dp.Basic_functions().array_fftx(E_c, self.time_window)
        freq_max = freq[np.where(abs(E_f) == max(abs(E_f)))][0]
        FWHM = dp.Basic_functions.array_FWHM(freq, abs(E_f))
        so_qnum = (freq_max-self.fast_freq_max)/FWHM
        return simps(E_c ** 2, dx=self.dx)-so_qnum/1000
        # t = self.time[(self.time>-1)&(self.time<1)]
        # x = self.get_combined_field(time_delays, t)
        # ts = abs(t[1]-t[0])
        # NFFT = int(2**(np.ceil(np.log2(len(t)))+2))
        # fs = 1/ts
        # x_normal = x-np.average(x)
        # f_stft, t_stft, zx = stft(x_normal, fs, nperseg = len(t)/2, 
        #                              noverlap=len(t)/2-1, window = ('gaussian', 12), nfft = NFFT)
        # # f, t, zx = self.fftx_stft(self.time, E)
        # # f_loc_mean = self.localf_so(self.time, f, zx)
        # f_ave = np.zeros_like(t)
        # # f_ave = np.average(f, axis = 0, weights = abs(zx))
        # for i in range(0, len(t)):
        #     # ampl = abs(zx)
        #     # ampl_sum = sum(ampl[:,i])
        #     f_ave[i] = np.average(f_stft, axis = 0, weights=abs(zx)[:,i])
        # f_ave = np.array(f_ave)
        # f_so_ave = np.mean(f_ave[(t >= -0.6) & (t <= 0.6)])
        
        # return 1/f_so_ave
    
    # def fftx_stft(self, t, x, pad = 2):
    #     ts = abs(t[1]-t[0])
    #     NFFT = int(2**(np.ceil(np.log2(len(t)))+pad))
    #     fs = 1/ts
    #     x_normal = x-np.average(x)
    #     freq_auto, t_auto, zx = stft(x_normal, fs, nperseg = len(t)/2, noverlap=len(t)/2-1, window = ('gaussian', 12), nfft = NFFT)
    #     return freq_auto, t_auto, zx
    
    # def localf_so(self, t, f, zx):
    #     f_ave = np.zeros_like(t)
    #     # f_ave = np.average(f, axis = 0, weights = abs(zx))
    #     for i in range(0, len(t)):
    #         # ampl = abs(zx)
    #         # ampl_sum = sum(ampl[:,i])
    #         f_ave[i] = np.average(f, axis = 0, weights=abs(zx)[:,i])
    #     f_ave = np.array(f_ave)
    #     f_so_ave = np.mean(f_ave[(t >= -0.6) & (t <= 0.6)])
    #     return f_so_ave
    
    # def locfreq2(self, t, f, zx):
        
    #     return np.sqrt(np.abs(
    #         np.gradient(np.gradient(E, t[1] - t[0]), t[1] - t[0]) / E
    #     ))
    # def jac_inegral_without_ampl_modulation(self,time_delays):
        
    #     E = self.get_combined_field(self,time_delays, self.time_window)
    
    #     return np.array([     
    #     -2. * simps(E * _.interp_field.derivative()(self.time_window - delay), dx=self.dx) 
    #     for delay, _ in zip(time_delays, self.pulses.values())
    # ])
    
    def local_minimization(self,initial_time_delays):
        result = minimize(
            self.inegral_without_ampl_modulation,
            initial_time_delays,
            #jac = jac_inegral_without_ampl_modulation,
            bounds=self.bounds,
            options={'maxiter': 1000},
        )
        # There is 2 decimal precision in the experiment 
        time_delays = np.round(result.x, 2)
        
        # new_delay = (0, time_delays[1], time_delays[2], time_delays[3])
        return self.inegral_without_ampl_modulation(time_delays), tuple(time_delays)

    def quant_so(self, time_delays):
        E_t = self.get_combined_field(time_delays, self.time)
        E_t_short = E_t[(self.time >= -0.6) & (self.time <= 0.6)]
        freq, E_f = dp.Basic_functions().array_fftx(E_t_short, self.time_window_raw)
        freq_max = freq[np.where(abs(E_f) == max(abs(E_f)))][0]
        FWHM = dp.Basic_functions.array_FWHM(freq, abs(E_f))
        so_qnum = (freq_max-self.fast_freq_max)/FWHM
        # d_qnum = abs(so_qnum-self.so_qnum)
        return so_qnum, time_delays
            
            
    def gradient_method(self):
        gradient_descent_results = set(self.local_minimization(_) for _ in tqdm(self.rand_initial_guess))
        gradient_descent_results = sorted(gradient_descent_results)
        intensity_without_ampl_modulation, all_time_delays = zip(*gradient_descent_results)
        return gradient_descent_results

    def get_err_combined_field(self,time_delays, time_window):
        return np.sqrt(sum(
            (_.err_upper_bound(time_window - delay) - _.err_lower_bound(time_window - delay) )** 2 
                for delay, _ in zip(time_delays, self.pulses.values())
            ))
   
    # def get_unique_filename(fname):
    #     return fname.format(datetime.now().strftime("_%m-%d_%H-%M"))    
    
    
class Calculate_optimal_phase:
    def __init__(self, path, fname0, fname1, filetype = '.dat'):
        self.pulses0 = {}
        self.pulses1 = {}
        self.max_ampl = []
        # make input variable into dict format
        
        #load data
        for name in fname0:
            # load all pulses for the given frequency
            self.filenames = glob.glob(path + name + '_*[0-9].dat')
            self.fnames_valid = []
            for j in range(0,len(self.filenames)):
                filenum = str(j+1)
                try:
                    data = np.loadtxt(path+'/'+name+'_'+filenum+filetype)
                except OSError:
                    continue
                else:
                    self.fnames_valid.append(path+'/'+name+'_'+filenum+filetype)
            self.all_data = [np.loadtxt(_).T[:-1] for _ in self.fnames_valid]
            self.times, self._individual_fields = zip(*self.all_data)
            self._individual_fields = np.array(self._individual_fields)
            time_mid = round((self.times[0][0]+self.times[0][-1])/2)
            time_diff = time_mid - 0
            self.time = self.times[0] - time_diff
            self.label = name
            self.field = np.mean(self._individual_fields, axis=0)
            # Substract the mean from experimental fields to compensate for parasitic offset
            self.field -= self.field.mean()
            # calculate the confidence interval
            confidence_level = 0.95
            self.sigma = stats.sem(self._individual_fields, axis=0)
            self.sigma = np.clip(self.sigma, 1e-20, self.sigma.max())
            self.err_lower_bound, self.err_upper_bound = stats.t.interval(confidence_level, 
                                                                          self._individual_fields.shape[0], 
                                                                          self.field, self.sigma,)
            abs_field = np.abs(self.field)
            self.max_ampl.append(abs_field.max())
            # Extract information for the combinatorial method     
            ampl_threshold = 0.4 * abs_field.max()
    
            indx = find_peaks(abs_field, height=ampl_threshold)[0]
            peaks_time = self.time[indx[1:-1]]
    
            self.half_period = Counter(np.diff(peaks_time)).most_common(1)[0][0]
        
                # Saving the data 
            self.pulses0[self.label] = Pulse(
                time = self.time, 
                field = self.field,
                _individual_fields = self._individual_fields,
                peaks_time = peaks_time,
                time_range = self.time[indx[1]:indx[-2]],
                half_period = self.half_period,
                # first save confidence interval as arrays
                err_lower_bound = self.err_lower_bound,
                err_upper_bound = self.err_upper_bound,
            )
        self.largest_freq_0 = self.label
            
        for name in fname1:
            # load all pulses for the given frequency
            self.filenames = glob.glob(path + name + '_*[0-9].dat')
            self.fnames_valid = []
            for j in range(0,len(self.filenames)):
                filenum = str(j+1)
                try:
                    data = np.loadtxt(path+'/'+name+'_'+filenum+filetype)
                except OSError:
                    continue
                else:
                    self.fnames_valid.append(path+'/'+name+'_'+filenum+filetype)
            self.all_data = [np.loadtxt(_).T[:-1] for _ in self.fnames_valid]
            self.times, self._individual_fields = zip(*self.all_data)
            self._individual_fields = np.array(self._individual_fields)
            time_mid = round((self.times[0][0]+self.times[0][-1])/2)
            time_diff = time_mid - 0
            self.time = self.times[0]-time_diff
            self.label = name
            self.field = np.mean(self._individual_fields, axis=0)
            # Substract the mean from experimental fields to compensate for parasitic offset
            self.field -= self.field.mean()
            # calculate the confidence interval
            confidence_level = 0.95
            self.sigma = stats.sem(self._individual_fields, axis=0)
            self.sigma = np.clip(self.sigma, 1e-20, self.sigma.max())
            self.err_lower_bound, self.err_upper_bound = stats.t.interval(confidence_level, 
                                                                          self._individual_fields.shape[0], 
                                                                          self.field, self.sigma,)
            abs_field = np.abs(self.field)
            self.max_ampl.append(abs_field.max())
            # Extract information for the combinatorial method     
            ampl_threshold = 0.4 * abs_field.max()
    
            indx = find_peaks(abs_field, height=ampl_threshold)[0]
            peaks_time = self.time[indx[1:-1]]
    
            self.half_period = Counter(np.diff(peaks_time)).most_common(1)[0][0]
        
                # Saving the data 
            self.pulses1[self.label] = Pulse(
                time = self.time, 
                field = self.field,
                _individual_fields = self._individual_fields,
                peaks_time = peaks_time,
                time_range = self.time[indx[1]:indx[-2]],
                half_period = self.half_period,
                # first save confidence interval as arrays
                err_lower_bound = self.err_lower_bound,
                err_upper_bound = self.err_upper_bound,
            )
        
        
        self.largest_freq_1 = self.label

        #  Checking whether the time axis coincide
        #assert all(np.allclose(time, data.time) for data in pulses.values()), \
            #"This workbook cannot be used since the time data is not syncronized"

        # saving time step
        self.dtime = self.time[1] - self.time[0]

        self.max_ampl = max(self.max_ampl)
            
            

        # Normalazing fields and interpolating
        for data0, data1 in zip(self.pulses0.values(), self.pulses1.values()):
            data0.field /= self.max_ampl
            data1.field /= self.max_ampl
            
            data0.interp_field = UnivariateSpline(
                data0.time, 
                data0.field,
                ext='zeros', 
                k=3, 
                s=0
            )
            data1.interp_field = UnivariateSpline(
                data1.time, 
                data1.field,
                ext='zeros', 
                k=3, 
                s=0
            )
            
        
            # interpolate confidence interval
            data0.err_lower_bound = UnivariateSpline(
                data0.time, 
                data0.err_lower_bound /self.max_ampl, 
                ext='zeros', 
                k=3, 
                s=0
            )
            
            data1.err_lower_bound = UnivariateSpline(
                data1.time, 
                data1.err_lower_bound /self.max_ampl, 
                ext='zeros', 
                k=3, 
                s=0
            )
            
            data0.err_upper_bound = UnivariateSpline(
                data0.time, 
                data0.err_upper_bound / self.max_ampl, 
                ext='zeros', 
                k=3, 
                s=0
            )
            
            data1.err_upper_bound = UnivariateSpline(
                data1.time, 
                data1.err_upper_bound / self.max_ampl, 
                ext='zeros', 
                k=3, 
                s=0
            )
         
            
            #print(data.field)
            #print(list(pulses))
        self.observational_window = -0.6, 0.6
        self.half_period = self.pulses0[self.largest_freq_0].half_period
        self.time_window, self.dx = np.linspace(-0.6, 0.6, 100, retstep=True)
        self.time_window_raw = self.time[(self.time >= -0.6) & (self.time <= 0.6)] 
        
        #generate initial guesses
        np.random.seed(3112022)

        self.bounds = [(_.peaks_time.min(), _.peaks_time.max()) for _ in self.pulses0.values()]
        
        #bounds_mid = [(_.peaks_time.min()+ _.peaks_time.max())/2 for _ in pulses.values()]
        #print(bounds_mid)
        self.rand_initial_guess_0 = np.array([
            np.random.uniform(_[0], _[1], 200 * cpu_count()) for _ in self.bounds
        ])
        self.rand_initial_guess_0 = set(tuple(_) for _ in self.rand_initial_guess_0.T)
        
        # self.rand_initial_guess_1 = np.array([
        #     np.random.uniform(_[0], _[1], 200 * cpu_count()) for _ in self.bounds
        # ])
        # self.rand_initial_guess_1 = set(tuple(_) for _ in self.rand_initial_guess_1.T)
        #print(rand_initial_guess)
        # Pool()
        pool = Pool(processes = int(cpu_count()/2))
      
        
        self.gradient_descent_results = set(pool.map(self.local_minimization, self.rand_initial_guess_0))
    
        
       
        self.gradient_descent_results = sorted(self.gradient_descent_results)
        self.intensity_without_ampl_modulation, self.all_time_delays = zip(*self.gradient_descent_results)
        # self.so_quant_results = []
        # for time_delays in self.all_time_delays:
        #     self.so_quant_results.append(self.quant_so(time_delays))
        # self.so_quant_results = sorted(self.so_quant_results, reverse = True)
        # self.so_qnums, self.all_time_delays = zip(*self.so_quant_results)
    
    def get_combined_field(self,time_delays, time_window):
        return sum(_.interp_field(time_window - delay) for delay, _ in zip(time_delays, self.pulses.values()))
    
    def get_combined_fields(self, time_delays, time_window):
        # time_delays_0 = time_delays[0]
        # time_delays_1 = time_delays[1]
        E0 = sum(_.interp_field(time_window - delay) for delay, _ in zip(time_delays, self.pulses0.values()))
        E1 = sum(_.interp_field(time_window - delay) for delay, _ in zip(time_delays, self.pulses1.values()))
        return E0, E1
   
    def quant_so(self, time_delays):
        E0_t, E1_t = self.get_combined_fields(time_delays, self.time)
        E0_t_short = E0_t[(self.time >= -0.6) & (self.time <= 0.6)]
        # E1_t_short = E1_t[(self.time >= -0.6) & (self.time <= 0.6)]
        freq0, E0_f = dp.Basic_functions().array_fftx(E0_t_short, self.time_window_raw)
        # freq1, E1_f = dp.Basic_functions().array_fftx(E1_t_short, self.time_window_raw)
        freq_max = freq0[np.where(abs(E0_f) == max(abs(E0_f)))][0]
        FWHM = dp.Basic_functions.array_FWHM(freq0, abs(E0_f))
        so_qnum = freq_max/FWHM
        # d_qnum = abs(so_qnum-self.so_qnum)
        return so_qnum, time_delays    
   
    def inegral_without_ampl_modulation(self,time_delays):
        return simps(self.get_combined_field(time_delays, self.time_window) ** 2, dx=self.dx)
    
    def get_normalized_discriminability(self, time_delays):
        E0, E1 = self.get_combined_fields(time_delays, self.time_window)
        J = simps(abs(E0-E1)**2)/simps(abs(E1)**2+abs(E0)**2)*2
        return 1/J
    
    def integral_inverse(self, time_delays):
        return 1/simps(self.get_combined_field(time_delays, self.time_window) ** 2, dx=self.dx)

    def jac_inegral_without_ampl_modulation(self,time_delays):
        
        E = self.get_combined_field(self,time_delays, self.time_window)
    
        return np.array([     
        -2. * simps(E * _.interp_field.derivative()(self.time_window - delay), dx=self.dx) 
        for delay, _ in zip(time_delays, self.pulses.values())
    ])
    def local_minimization(self,initial_time_delays):
        # initial_time_delays_0 = initial_time_delays[0]
        # initial_time_delays_0 = np.array(initial_time_delays_0)
        # initial_time_delays_1 = initial_time_delays[1]
        # initial_time_delays_1 = np.array(initial_time_delays_1)
        initial_time_delays = np.array(initial_time_delays)
        result = minimize(
            self.get_normalized_discriminability,
            initial_time_delays,
            #jac = jac_inegral_without_ampl_modulation,
            bounds=self.bounds,
            options={'maxiter': 1000},
        )
        # There is 2 decimal precision in the experiment 
        time_delays = np.round(result.x, 2)
        return self.get_normalized_discriminability(time_delays), tuple(time_delays)

    def gradient_method(self):
        gradient_descent_results = set(self.local_minimization(_) for _ in tqdm(self.rand_initial_guess))
        gradient_descent_results = sorted(gradient_descent_results)
        intensity_without_ampl_modulation, all_time_delays = zip(*gradient_descent_results)
        return gradient_descent_results

    def get_err_combined_field(self,time_delays, time_window):
        return np.sqrt(sum(
            (_.err_upper_bound(time_window - delay) - _.err_lower_bound(time_window - delay) )** 2 
                for delay, _ in zip(time_delays, self.pulses0.values())
            ))
   
    # def get_unique_filename(fname):
    #     return fname.format(datetime.now().strftime("_%m-%d_%H-%M"))    
    
class Calculate_optimal_phase_varwindow:
    def __init__(self, path, fname0, fname1, window, filetype = '.dat'):
        self.pulses0 = {}
        self.pulses1 = {}
        self.max_ampl = []
        self.designate_window = window
        # make input variable into dict format
        
        #load data
        for name in fname0:
            # load all pulses for the given frequency
            self.filenames = glob.glob(path + name + '_*[0-9].dat')
            self.fnames_valid = []
            for j in range(0,len(self.filenames)):
                filenum = str(j+1)
                try:
                    data = np.loadtxt(path+'/'+name+'_'+filenum+filetype)
                except OSError:
                    continue
                else:
                    self.fnames_valid.append(path+'/'+name+'_'+filenum+filetype)
            self.all_data = [np.loadtxt(_).T[:-1] for _ in self.fnames_valid]
            self.times, self._individual_fields = zip(*self.all_data)
            self._individual_fields = np.array(self._individual_fields)
            time_mid = round((self.times[0][0]+self.times[0][-1])/2)
            time_diff = time_mid - 0
            self.time = self.times[0] - time_diff
            self.label = name
            self.field = np.mean(self._individual_fields, axis=0)
            # Substract the mean from experimental fields to compensate for parasitic offset
            self.field -= self.field.mean()
            # calculate the confidence interval
            confidence_level = 0.95
            self.sigma = stats.sem(self._individual_fields, axis=0)
            self.sigma = np.clip(self.sigma, 1e-20, self.sigma.max())
            self.err_lower_bound, self.err_upper_bound = stats.t.interval(confidence_level, 
                                                                          self._individual_fields.shape[0], 
                                                                          self.field, self.sigma,)
            abs_field = np.abs(self.field)
            self.max_ampl.append(abs_field.max())
            # Extract information for the combinatorial method     
            ampl_threshold = 0.4 * abs_field.max()
    
            indx = find_peaks(abs_field, height=ampl_threshold)[0]
            peaks_time = self.time[indx[1:-1]]
    
            self.half_period = Counter(np.diff(peaks_time)).most_common(1)[0][0]
        
                # Saving the data 
            self.pulses0[self.label] = Pulse(
                time = self.time, 
                field = self.field,
                _individual_fields = self._individual_fields,
                peaks_time = peaks_time,
                time_range = self.time[indx[1]:indx[-2]],
                half_period = self.half_period,
                # first save confidence interval as arrays
                err_lower_bound = self.err_lower_bound,
                err_upper_bound = self.err_upper_bound,
            )
        self.largest_freq_0 = self.label
            
        for name in fname1:
            # load all pulses for the given frequency
            self.filenames = glob.glob(path + name + '_*[0-9].dat')
            self.fnames_valid = []
            for j in range(0,len(self.filenames)):
                filenum = str(j+1)
                try:
                    data = np.loadtxt(path+'/'+name+'_'+filenum+filetype)
                except OSError:
                    continue
                else:
                    self.fnames_valid.append(path+'/'+name+'_'+filenum+filetype)
            self.all_data = [np.loadtxt(_).T[:-1] for _ in self.fnames_valid]
            self.times, self._individual_fields = zip(*self.all_data)
            self._individual_fields = np.array(self._individual_fields)
            time_mid = round((self.times[0][0]+self.times[0][-1])/2)
            time_diff = time_mid - 0
            self.time = self.times[0] - time_diff
            # time_mid = self.times[0][round((self.times[0][0]+self.times[0][-1])/2)]
            # time_diff = time_mid - 0
            # self.time = self.times[0]
            self.label = name
            self.field = np.mean(self._individual_fields, axis=0)
            # Substract the mean from experimental fields to compensate for parasitic offset
            self.field -= self.field.mean()
            # calculate the confidence interval
            confidence_level = 0.95
            self.sigma = stats.sem(self._individual_fields, axis=0)
            self.sigma = np.clip(self.sigma, 1e-20, self.sigma.max())
            self.err_lower_bound, self.err_upper_bound = stats.t.interval(confidence_level, 
                                                                          self._individual_fields.shape[0], 
                                                                          self.field, self.sigma,)
            abs_field = np.abs(self.field)
            self.max_ampl.append(abs_field.max())
            # Extract information for the combinatorial method     
            ampl_threshold = 0.4 * abs_field.max()
    
            indx = find_peaks(abs_field, height=ampl_threshold)[0]
            peaks_time = self.time[indx[1:-1]]
    
            self.half_period = Counter(np.diff(peaks_time)).most_common(1)[0][0]
        
                # Saving the data 
            self.pulses1[self.label] = Pulse(
                time = self.time, 
                field = self.field,
                _individual_fields = self._individual_fields,
                peaks_time = peaks_time,
                time_range = self.time[indx[1]:indx[-2]],
                half_period = self.half_period,
                # first save confidence interval as arrays
                err_lower_bound = self.err_lower_bound,
                err_upper_bound = self.err_upper_bound,
            )
        
        
        self.largest_freq_1 = self.label

        #  Checking whether the time axis coincide
        #assert all(np.allclose(time, data.time) for data in pulses.values()), \
            #"This workbook cannot be used since the time data is not syncronized"

        # saving time step
        self.dtime = self.time[1] - self.time[0]

        self.max_ampl = max(self.max_ampl)
            
            

        # Normalazing fields and interpolating
        for data0, data1 in zip(self.pulses0.values(), self.pulses1.values()):
            data0.field /= self.max_ampl
            data1.field /= self.max_ampl
            
            data0.interp_field = UnivariateSpline(
                data0.time, 
                data0.field,
                ext='zeros', 
                k=3, 
                s=0
            )
            data1.interp_field = UnivariateSpline(
                data1.time, 
                data1.field,
                ext='zeros', 
                k=3, 
                s=0
            )
            
        
            # interpolate confidence interval
            data0.err_lower_bound = UnivariateSpline(
                data0.time, 
                data0.err_lower_bound /self.max_ampl, 
                ext='zeros', 
                k=3, 
                s=0
            )
            
            data1.err_lower_bound = UnivariateSpline(
                data1.time, 
                data1.err_lower_bound /self.max_ampl, 
                ext='zeros', 
                k=3, 
                s=0
            )
            
            data0.err_upper_bound = UnivariateSpline(
                data0.time, 
                data0.err_upper_bound / self.max_ampl, 
                ext='zeros', 
                k=3, 
                s=0
            )
            
            data1.err_upper_bound = UnivariateSpline(
                data1.time, 
                data1.err_upper_bound / self.max_ampl, 
                ext='zeros', 
                k=3, 
                s=0
            )
         
            
            #print(data.field)
            #print(list(pulses))
        # self.observational_window = -0.6, 0.6
        # self.half_period = self.pulses0[self.largest_freq_0].half_period
        self.time_window, self.dx = np.linspace(self.designate_window[0], self.designate_window[1], 100, retstep=True)
        self.time_window_raw = self.time[(self.time >= self.designate_window[0]) & (self.time <= self.designate_window[1])] 
        
        #generate initial guesses
        np.random.seed(3112022)

        self.bounds = [(_.peaks_time.min(), _.peaks_time.max()) for _ in self.pulses0.values()]
        
        #bounds_mid = [(_.peaks_time.min()+ _.peaks_time.max())/2 for _ in pulses.values()]
        #print(bounds_mid)
        self.rand_initial_guess_0 = np.array([
            np.random.uniform(_[0], _[1], 200 * cpu_count()) for _ in self.bounds
        ])
        self.rand_initial_guess_0 = self.rand_initial_guess_0.T
        for i in range(0,len(self.rand_initial_guess_0)):
            self.rand_initial_guess_0[i,0] = 0
            
        self.rand_initial_guess_0 = set(tuple(_) for _ in self.rand_initial_guess_0)
        
            
       
        pool = Pool(processes = int(cpu_count()))
        
        self.gradient_descent_results = set(pool.map(self.local_minimization, self.rand_initial_guess_0))
        # self.gradient_descent_results = set(self.local_minimization(_) for _ in (self.rand_initial_guess_0))
    
        self.gradient_descent_results = sorted(self.gradient_descent_results)
        self.intensity_without_ampl_modulation, self.all_time_delays = zip(*self.gradient_descent_results)
        self.best_J = 1/self.get_normalized_discriminability(self.all_time_delays[0])
        self.localf_E0, self.localf_E1 = self.local_f_so(self.all_time_delays[0], self.time)
        self.localf_ave = (np.mean(self.localf_E0)+np.mean(self.localf_E1))/2
        self.field_E0, self.field_E1 = self.get_combined_fields(self.all_time_delays[0], self.time)
        
    
    def get_combined_field(self,time_delays, time_window):
        return sum(_.interp_field(time_window - delay) for delay, _ in zip(time_delays, self.pulses.values()))
    
    def get_combined_fields(self, time_delays, time_window):
        # time_delays_0 = time_delays[0]
        # time_delays_1 = time_delays[1]
        # new_delay = (0, time_delays[1], time_delays[2], time_delays[3])
        # E0 = sum(_.interp_field(time_window - delay) for delay, _ in zip(new_delay, self.pulses0.values()))
        # E1 = sum(_.interp_field(time_window - delay) for delay, _ in zip(new_delay, self.pulses1.values()))
        E0 = sum(_.interp_field(time_window - delay) for delay, _ in zip(time_delays, self.pulses0.values()))
        E1 = sum(_.interp_field(time_window - delay) for delay, _ in zip(time_delays, self.pulses1.values()))
        return E0, E1
    
    def inegral_without_ampl_modulation(self,time_delays):
        return simps(self.get_combined_field(time_delays, self.time_window) ** 2, dx=self.dx)
    
    def get_normalized_discriminability(self, time_delays):
        E0, E1 = self.get_combined_fields(time_delays, self.time_window)
        J = simps(abs(E0-E1)**2)/simps(abs(E1)**2+abs(E0)**2)*2
        # J = simps(abs(E0-E1)**2)
        return 1/J
    
    def local_f_so(self, time_delays, time_window):
        E0 = sum(_.interp_field(time_window - delay) for delay, _ in zip(time_delays, self.pulses0.values()))
        E1 = sum(_.interp_field(time_window - delay) for delay, _ in zip(time_delays, self.pulses1.values()))
        t_stft, localf_E0 = fdp.localf_stft(time_window, E0)
        t_stft, localf_E1 = fdp.localf_stft(time_window, E1)
        localf_E0_so = localf_E0[(t_stft >= self.designate_window[0]) & (t_stft <= self.designate_window[1])]
        localf_E1_so = localf_E1[(t_stft >= self.designate_window[0]) & (t_stft <= self.designate_window[1])]
        return localf_E0_so, localf_E1_so
        
        
    def jac_inegral_without_ampl_modulation(self,time_delays):
        
        E = self.get_combined_field(self,time_delays, self.time_window)
    
        return np.array([     
        -2. * simps(E * _.interp_field.derivative()(self.time_window - delay), dx=self.dx) 
        for delay, _ in zip(time_delays, self.pulses.values())
    ])
    def local_minimization(self,initial_time_delays):
        # initial_time_delays_0 = initial_time_delays[0]
        # initial_time_delays_0 = np.array(initial_time_delays_0)
        # initial_time_delays_1 = initial_time_delays[1]
        # initial_time_delays_1 = np.array(initial_time_delays_1)
        initial_time_delays = np.array(initial_time_delays)
        result = minimize(
            self.get_normalized_discriminability,
            initial_time_delays,
            #jac = jac_inegral_without_ampl_modulation,
            bounds=self.bounds,
            options={'maxiter': 1000},
        )
        # There is 2 decimal precision in the experiment 
        time_delays = np.round(result.x, 2)
        # new_delays = (0, time_delays[1], time_delays[2], time_delays[3])
        
        return self.get_normalized_discriminability(time_delays), tuple(time_delays)

    def gradient_method(self):
        gradient_descent_results = set(self.local_minimization(_) for _ in tqdm(self.rand_initial_guess))
        gradient_descent_results = sorted(gradient_descent_results)
        intensity_without_ampl_modulation, all_time_delays = zip(*gradient_descent_results)
        return gradient_descent_results

    def get_err_combined_field(self,time_delays, time_window):
        return np.sqrt(sum(
            (_.err_upper_bound(time_window - delay) - _.err_lower_bound(time_window - delay) )** 2 
                for delay, _ in zip(time_delays, self.pulses0.values())
            ))
   
    # def get_unique_filename(fname):
    #     return fname.format(datetime.now().strftime("_%m-%d_%H-%M")) 
        
if __name__ == '__main__':
    J_max = []
    path = 'C:/Users/Admin/Dropbox/SO project/very different fields/'
    fname0 = ['H_N0S_F50', 'H_N0S_F60', 'H_N0S_F70', 'H_N0S_F80']
    fname1 = ['H_N2S_F50', 'H_N2S_F60', 'H_N2S_F70', 'H_N2S_F80']
    # fname0 = ['CX0_F05', 'CX0_F06', 'CX0_F07', 'CX0_F08']
    # fname1 = ['CY0_F05', 'CY0_F06', 'CY0_F07', 'CY0_F08']
    half_window_len = np.arange(0.2, 4, step = 0.1)
    localf_ave = []
    localf_E0 = []
    localf_E1 = []
    plt.figure()
    for idx, h_win in enumerate(half_window_len):
        so_optJ_go = Calculate_optimal_phase_varwindow(path, fname0, fname1, 
                                                        window = (-h_win, h_win))
        J_max.append(so_optJ_go.intensity_without_ampl_modulation[0])
        localf_ave.append(so_optJ_go.localf_ave)
        localf_E0.append(so_optJ_go.localf_E0)
        localf_E1.append(so_optJ_go.localf_E1)       
        plt.plot(so_optJ_go.time[(so_optJ_go.time>=-5)&(so_optJ_go.time<=5)], 
                 so_optJ_go.field_E0[(so_optJ_go.time>=-5)&(so_optJ_go.time<=5)])
        
        print('window #'+str(idx)+'calculated')
    J_max = 1/np.array(J_max)
    localf_ave = np.array(localf_ave)
    zero_filler = np.zeros_like(J_max)
    spec = np.stack((2*half_window_len, np.log10(J_max), zero_filler), axis = 1)
    np.savetxt(path+'Jmax_1.dat', spec)
    plt.figure()
    plt.plot(2*half_window_len, np.log10(J_max))
    
    plt.xlabel('time window length (ps)')
    plt.ylabel('$\log_{10}(J)$')
    plt.figure()
    plt.plot(2*half_window_len, localf_ave)
    plt.xlabel('time window length (ps)')
    plt.ylabel('average local frequency (THz)')
