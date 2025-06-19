# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 21:13:40 2023

@author: ppsap
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:49:14 2022

@author: ppsapphire
"""

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import tkinter as tk
import dataprocessor as dp
import glob
import plottool_v6 as pt
import copy

from scipy import stats
from scipy.signal import find_peaks, hilbert, windows
from scipy.integrate import simps
from scipy.interpolate import UnivariateSpline, InterpolatedUnivariateSpline
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize
from scipy.signal import stft, spectrogram


from collections import Counter
from itertools import product
from dataclasses import dataclass
from threading import Thread

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk
from tkinter import ttk

from multiprocessing import Pool, cpu_count, Process, Queue
from tqdm.notebook import tqdm
import configparser as cfp

import SO_calculator as socal




from datetime import datetime




class SO_init:
    '''generate a GUI for superoscillation related applications'''
    
    def __init__(self, root_window, datapath, x_avg, t_avg):
        '''
        

        Parameters
        ----------
        root_window : Tk Type
            DESCRIPTION. main window for this application
        x_avg : Dict TYPE
            DESCRIPTION. Dictionary that contains averaged x values to be used in calculation
        t_avg : Dict TYPE
            DESCRIPTION. Dictionary that contains averaged time values to be used in calculation

        Returns
        -------
        None.

        '''
        
        
        '''global varibales'''
        self.root = root_window
        # store the values of x and t in deepcopy form so any modification to x and t dicts won't influence original dicts
        self.x = copy.deepcopy(x_avg)
        self.t = copy.deepcopy(t_avg)
        # also get a backup copy in case we want to restore to original values after certain manipulations
        self.x_backup = copy.deepcopy(x_avg)
        self.t_backup = copy.deepcopy(t_avg)
        
        
        
        self.path = datapath + '/'
        self.filenames = sorted(list(t_avg)) # get the dict keys of x and t 
        self.status = tk.StringVar(value = 'ready') #s tatus bar textvaribale for displaying status of the application
        self.listbox_filenames_content = tk.StringVar(value = self.filenames) 
        self.J = {}
        
        self.J_window = []
        self.J_conv = {}
        self.J_base = {}
        self.J_log = {}
        self.tw_len = {}
        self.f_J = {}
        self.n_dis = 0
        self.D_J = tk.StringVar(value = 'none')
        
        # textvaribales for displaying names of selected data sets
        self.selection_wf_0 = tk.StringVar(value = 'not selected')
        self.selection_wf_1 = tk.StringVar(value = 'not selected')
        self.selection_so_0 = tk.StringVar(value = 'not selected')
        self.selection_so_1 = tk.StringVar(value = 'not selected')
        
        # text variables for storing and dispaly time delays
        self.td1 = tk.StringVar()
        self.td2 = tk.StringVar()
        self.td3 = tk.StringVar()
        
        # textvariables for storing and display newly calculated phases
        self.ph1 = tk.StringVar()
        self.ph2 = tk.StringVar()
        self.ph3 = tk.StringVar()
        
        # shortcut handle to use Basic_functions module
        self.bf = dp.Basic_functions()
        
        # textvariable to store and display the starting index of superoscillation results to display
        self.n0 = tk.IntVar(value=0)
        
        # variables to store index of selected data sets
        # self.idx_so_0 = (0,)
        # self.idx_so_1 = (0,)
        # self.idx_wf_0 = (0,)
        # self.idx_wf_1 = (0,)
        self.idx_so_0 = ()
        self.idx_so_1 = ()
        self.idx_wf_0 = ()
        self.idx_wf_1 = ()
        
        
        # try to load pre-configurated values to textvariables that are often changed, this is to save time of entering these values 
        # when the window restarts
        try:
            open('so_init.ini')
        except:
            self.timewindow_bot = tk.DoubleVar()
            self.timewindow_top = tk.DoubleVar()
            
            self.phase_05 = tk.DoubleVar(value=-10)
            self.phase_06 = tk.DoubleVar(value=-133)
            self.phase_07 = tk.DoubleVar(value=-114)
            self.phase_08 = tk.DoubleVar(value=-26)
            
            
        else:
            self.timewindow_bot = tk.DoubleVar()
            self.timewindow_top = tk.DoubleVar()
            
            self.phase_05 = tk.DoubleVar()
            self.phase_06 = tk.DoubleVar()
            self.phase_07 = tk.DoubleVar()
            self.phase_08 = tk.DoubleVar()
        
            so_config = cfp.ConfigParser().read('so_tools.ini')
            so_config_keys = list(so_config['so_init'])
            
            if 'timewindow_bot' in so_config_keys:
                self.timewindow_bot.set(float(so_config['so_init']['timewindow_bot']))
            if 'timewindow_top' in so_config_keys:
                self.timewindow_top.set(float(so_config['so_init']['timewindow_top']))
            if 'phase_05' in so_config_keys:
                self.phase_05.set(float(so_config['so_init']['phase_05']))
            if 'phase_06' in so_config_keys:
                self.phase_05.set(float(so_config['so_init']['phase_06']))
            if 'phase_07' in so_config_keys:
                self.phase_05.set(float(so_config['so_init']['phase_07']))
            if 'phase_08' in so_config_keys:
                self.phase_05.set(float(so_config['so_init']['phase_08']))
            
            
        '''main windows'''
        
        self.window_so_calculate = ttk.Labelframe(self.root, text = 'Calculate SO phases')
        self.window_so_discriminate = ttk.Labelframe(self.root, text = 'Calculate discrimination')
        self.window_so_plot = ttk.Labelframe(self.root, text = 'Plot results')
        self.window_so_status = ttk.Labelframe(self.root, text = 'status: ')
        self.window_so_calculate.grid(column = 0, row = 0)
        self.window_so_discriminate.grid(column = 0, row = 1)
        self.window_so_plot.grid(column = 1, row = 0, columnspan=2, rowspan = 2)
        self.window_so_status.grid(column = 0, row = 3, sticky = 'w')
        
        '''window content for window_so_calcualtet'''
        self.listbox_filenames = tk.Listbox(self.window_so_calculate, listvariable = self.listbox_filenames_content, 
                                            selectmode = 'extended', width = 40)
        label_wf_0 = ttk.Label(self.window_so_calculate, text = 'waveforms_0: \n')
        label_wf_0_content = ttk.Label(self.window_so_calculate, textvariable = self.selection_wf_0)
        label_wf_1 = ttk.Label(self.window_so_calculate, text = 'waveforms_1: \n')
        label_wf_1_content = ttk.Label(self.window_so_calculate, textvariable = self.selection_wf_1)
        label_so_0 = ttk.Label(self.window_so_calculate, text = 'superosc_0: \n')
        label_so_0_content = ttk.Label(self.window_so_calculate, textvariable = self.selection_so_0)
        label_so_1 = ttk.Label(self.window_so_calculate, text = 'superosc_1: \n')
        label_so_1_content = ttk.Label(self.window_so_calculate, textvariable = self.selection_so_1)
        self.button_select_wf_0 = ttk.Button(self.window_so_calculate, text = 'select WF 0', command = self.select_wf_0)
        self.button_select_wf_1 = ttk.Button(self.window_so_calculate, text = 'select WF 1', command = self.select_wf_1)
        self.button_select_so_0 = ttk.Button(self.window_so_calculate, text = 'select SO 0', command = self.select_so_0)
        self.button_select_so_1 = ttk.Button(self.window_so_calculate, text = 'select SO 1', command = self.select_so_1)
        button_combine_wf = ttk.Button(self.window_so_calculate, text = 'combine wf 1 & 2', command = self.combine_waveforms_12)
        self.button_calculate_so = ttk.Button(self.window_so_calculate, text = 'Calculate SO', command = self.start_thread_cal_so)
        # button_calculate_so = ttk.Button(self.window_so_calculate, text = 'Calculate SO', command = self.calculate_so_phases)
        button_local_f = ttk.Button(self.window_so_calculate, text = 'Local freq', command = self.Local_f)
        button_calculate_phase = ttk.Button(self.window_so_calculate, text = 'show new phases', command = self.calculate_phase)
        
        entry_lb = ttk.Entry(self.window_so_calculate,textvariable=self.n0)
        # entry_up = ttk.Entry(self.window,textvariable=self.n1)
        label_lb = ttk.Label(self.window_so_calculate,text = 'Pick the 1st SO spectrum to display: ')
        button_show_timedelays = ttk.Button(self.window_so_calculate,text='show time delays', command = self.display_timedelays)
        label_timedelay1 = ttk.Label(self.window_so_calculate,text = 'time delay 1: ')
        label_timedelay2 = ttk.Label(self.window_so_calculate,text = 'time delay 2: ')
        label_timedelay3 = ttk.Label(self.window_so_calculate,text = 'time delay 3: ')
        label_timedelay1_content = ttk.Label(self.window_so_calculate,textvariable = self.td1)
        label_timedelay2_content = ttk.Label(self.window_so_calculate,textvariable = self.td2)
        label_timedelay3_content = ttk.Label(self.window_so_calculate,textvariable = self.td3)
        entry_05_phase = ttk.Entry(self.window_so_calculate,textvariable=self.phase_05)
        entry_06_phase = ttk.Entry(self.window_so_calculate,textvariable=self.phase_06)
        entry_07_phase = ttk.Entry(self.window_so_calculate,textvariable=self.phase_07)
        entry_08_phase = ttk.Entry(self.window_so_calculate,textvariable=self.phase_08)
        label_ph1 = ttk.Label(self.window_so_calculate,textvariable=self.ph1)
        label_ph2 = ttk.Label(self.window_so_calculate,textvariable=self.ph2)
        label_ph3 = ttk.Label(self.window_so_calculate,textvariable=self.ph3)
        label_05 = ttk.Label(self.window_so_calculate,text='0.5: ')
        label_06 = ttk.Label(self.window_so_calculate,text='0.6: ')
        label_07 = ttk.Label(self.window_so_calculate,text='0.7: ')
        label_08 = ttk.Label(self.window_so_calculate,text='0.8: ')
        label_cal_phase = ttk.Label(self.window_so_calculate,text='output phases: ')
        
        button_plot_so = ttk.Button(self.window_so_calculate, text = 'plot results', command = self.plot_predictions)
        
        button_sim_so = ttk.Button(self.window_so_calculate, text = 'simulate SO', command = self.simulate_so)
        
        button_optimal = ttk.Button(self.window_so_calculate, text = 'calculate max ND', command = self.start_thread_cal_optJ)
        button_compare_OP_SO = ttk.Button(self.window_so_calculate, text = 'compare SO with max ND', command = self.compare_so_op)
        button_shift_wfs = ttk.Button(self.window_so_calculate, text = 'shift waveforms', command = self.shift_wfs)
        
        button_save_config = ttk.Button(self.window_so_calculate, text = 'make default inputs', command = self.build_so_config)
        button_plotall = ttk.Button(self.window_so_calculate, text = 'plot all', command = self.plot_all)
        button_c_so = ttk.Button(self.window_so_calculate, text = 'auto SO', command = self.find_closest_so)
        
        
        self.listbox_filenames.grid(column = 0, row = 0, rowspan = 4, columnspan = 4)
        self.button_select_wf_0.grid(column = 5, row = 0)
        self.button_select_wf_1.grid(column = 5, row = 1)
        self.button_select_so_0.grid(column = 5, row = 2)
        self.button_select_so_1.grid(column = 5, row = 3)
        
        
        
        
        label_wf_0.grid(column = 0, row = 4)
        label_wf_0_content.grid(column = 0, row = 5)
        label_wf_1.grid(column = 1, row = 4)
        label_wf_1_content.grid(column = 1, row = 5)
        label_so_0.grid(column = 2, row = 4)
        label_so_0_content.grid(column = 2, row = 5)
        label_so_1.grid(column = 3, row = 4)
        label_so_1_content.grid(column = 3, row = 5)
        
        button_combine_wf.grid(column = 0, row = 6)
        self.button_calculate_so.grid(column = 1, row = 6)
        button_local_f.grid(column = 2, row = 6)
        button_show_timedelays.grid(column=3,row = 6)
        button_calculate_phase.grid(column = 4, row = 6)
        button_plot_so.grid(column = 5, row = 6)
        
        label_lb.grid(column=0,row=7)
        entry_lb.grid(column=1,row=7)
        label_timedelay1.grid(column=0,row=8)
        label_timedelay1_content.grid(column=1,row=8)
        label_timedelay2.grid(column=0,row=8)
        label_timedelay2_content.grid(column=1,row=8)
        label_timedelay3.grid(column=0,row=8)
        label_timedelay3_content.grid(column=1,row=8)
        label_05.grid(column=0,row=9)
        entry_05_phase.grid(column=1,row=9)
        label_06.grid(column=2,row=9)
        entry_06_phase.grid(column=3,row=9)
        label_07.grid(column=0,row=10)
        entry_07_phase.grid(column=1,row=10)
        label_08.grid(column=2,row=10)
        entry_08_phase.grid(column=3,row=10)
        label_cal_phase.grid(column=0,row=11)
        label_ph1.grid(column=1,row=11)
        label_ph2.grid(column=1,row=12)
        label_ph3.grid(column=1,row=13)
        
        button_shift_wfs.grid(column = 0, row = 14, sticky = 'w')
        button_sim_so.grid(column = 1, row = 14, sticky = 'w')
        button_optimal.grid(column = 2, row = 14, sticky = 'w')
        button_compare_OP_SO.grid(column = 3, row = 14, sticky = 'w')
        button_save_config.grid(column = 0, row = 15, sticky = 'w')
        button_plotall.grid(column = 1, row = 15, sticky = 'w')
        button_c_so.grid(column = 2, row = 15, sticky = 'w')
        
        
        '''window content for window_so_discriminate'''
        
        label_tw_bot = ttk.Label(self.window_so_discriminate, text = 'Enter lower limit of timewindow: ')
        label_tw_top = ttk.Label(self.window_so_discriminate, text = 'Enter upper limit of timewindow: ')
        entry_timewindow_bot = ttk.Entry(self.window_so_discriminate, textvariable = self.timewindow_bot)
        entry_timewindow_top = ttk.Entry(self.window_so_discriminate, textvariable = self.timewindow_top)
        button_cal_discrim = ttk.Button(self.window_so_discriminate, text = 'calculate ND', command = self.discrimin_normal)        
        button_D_discrim = ttk.Button(self.window_so_discriminate, text = 'calculate D ND', command = self.D_discrimin)
        button_fd_discrim = ttk.Button(self.window_so_discriminate, text = 'calculate fd ND', command = self.discrimin_normal_fd)
        button_conv_discrim = ttk.Button(self.window_so_discriminate, text = 'calculate conv D', command = self.discrimin_nd_shift)
        button_J_all = ttk.Button(self.window_so_discriminate, text = 'plot J all delay', command = self.J_all_optimal)
        button_J_varwindow = ttk.Button(self.window_so_discriminate, text = 'cal J varwindow', command = self.start_thread_maxJ_var_window)
        label_D_discrim = ttk.Label(self.window_so_discriminate, text = 'Advantage ratio: ')
        label_D_discrim_value = ttk.Label(self.window_so_discriminate, textvariable = self.D_J)
        self.discrimin_type = tk.StringVar(value = 'log scale')
        combobox_discrm_type = ttk.Combobox(self.window_so_discriminate, textvariable = self.discrimin_type)
        combobox_discrm_type['values'] = ['normalized', 'original', 'log scale']
        
        
        
        label_tw_bot.grid(column = 0, row = 0, sticky = 'w')
        entry_timewindow_bot.grid(column = 1, row = 0, sticky = 'w')
        label_tw_top.grid(column = 2, row = 0, sticky = 'w')
        entry_timewindow_top.grid(column = 3, row = 0, sticky = 'w')
        combobox_discrm_type.grid(column = 0, row = 1, sticky = 'w')
        button_cal_discrim.grid(column = 1, row = 1, sticky = 'w')
        button_D_discrim.grid(column = 2, row = 1, sticky = 'w')
        button_fd_discrim.grid(column = 3, row = 1, sticky = 'w')
        button_conv_discrim.grid(column = 4, row = 1, sticky = 'w')
        button_J_varwindow.grid(column = 0, row = 2, sticky = 'w')
        label_D_discrim.grid(column = 0, row = 3, sticky = 'w')
        label_D_discrim_value.grid(column = 1, row = 3, sticky = 'w')
        button_J_all.grid(column = 0, row = 4, sticky = 'w')
        
        
        '''window content for window_so_plot'''
        
        fig0 = Figure(figsize=(9,6), dpi = 100)
        self.canvas0 = FigureCanvasTkAgg(figure=fig0,master=self.window_so_plot)
        self.canvas0.get_tk_widget().grid(column=0,row=0,columnspan=6,rowspan=6,sticky='nsew')
        self.canvas0.figure, self.ax = plt.subplots(2,1,figsize=(4,3))
        plt.tight_layout()
        self.canvas0.draw()
        
        '''window content for window_so_status'''
        self.label_status = ttk.Label(self.window_so_status, textvariable = self.status, font=('Times', 15))
        self.progressbar = ttk.Progressbar(self.window_so_status, orient='horizontal', length=50, mode='indeterminate')
        # self.progressbar['maximum'] = 100
        self.label_status.grid(column = 0, row = 0, columnspan = 3)
        self.progressbar.grid(column = 4, row = 0, columnspan = 6)
        
        
    
    '''functions for button activities''' 
    
    
    def build_so_config(self):
        so_config = cfp.ConfigParser()
        so_config['so_init'] = {'timewindow_bot': str(self.timewindow_bot.get()),
                                'timewindow_top': str(self.timewindow_top.get()),
                                'phase_05': str(self.phase_05.get()), 
                                'phase_06': str(self.phase_06.get()),
                                'phase_07': str(self.phase_07.get()),
                                'phase_08': str(self.phase_08.get())
                                }
        with open('so_init.ini','w') as so_init:
            so_config.write(so_init)
    def select_wf_0(self):
        self.filenames_wf_0 = []
        self.interpo_wf0 = {}
        self.idx_wf_0 = self.listbox_filenames.curselection()
        selection_fields = ''
        for i in self.idx_wf_0:
            selection_fields = selection_fields + self.filenames[i] + '\n'
            self.filenames_wf_0.append(self.filenames[i])
        self.selection_wf_0.set(selection_fields)
        for name in self.filenames_wf_0:
            self.interpo_wf0[name] = UnivariateSpline(self.t[name], self.x[name], ext = 'zeros', k = 3, s = 0)
            
    def select_wf_1(self):
        self.filenames_wf_1 = []
        self.interpo_wf1 = {}
        self.idx_wf_1 = self.listbox_filenames.curselection()
        selection_fields = ''
        for i in self.idx_wf_1:
            selection_fields = selection_fields + self.filenames[i] + '\n'
            self.filenames_wf_1.append(self.filenames[i])
        self.selection_wf_1.set(selection_fields)
        for name in self.filenames_wf_1:
            self.interpo_wf1[name] = UnivariateSpline(self.t[name], self.x[name], ext = 'zeros', k = 3, s = 0)
    
    def select_so_0(self):
        self.filenames_so_0 = []
        self.idx_so_0 = self.listbox_filenames.curselection()
        selection_fields = ''
        for i in self.idx_so_0:
            selection_fields = selection_fields + self.filenames[i] + '\n'
            self.filenames_so_0.append(self.filenames[i])
        self.selection_so_0.set(selection_fields)
        
    def select_so_1(self):
        self.filenames_so_1 = []
        self.idx_so_1 = self.listbox_filenames.curselection()
        selection_fields = ''
        for i in self.idx_so_1:
            selection_fields = selection_fields + self.filenames[i] + '\n'
            self.filenames_so_1.append(self.filenames[i])
        self.selection_so_1.set(selection_fields)

            
        
    def combine_waveforms_12(self):
        self.combined_wf = {}
        if len(self.idx_wf_0) == len(self.idx_wf_1):
            for i in range(0,len(self.idx_wf_0)):
                self.combined_wf[self.filenames[self.idx_wf_0[i]]] = np.stack((self.t[self.filenames[self.idx_wf_0[i]]], 
                                                                                              self.x[self.filenames[self.idx_wf_0[i]]], 
                                                                                              self.x[self.filenames[self.idx_wf_1[i]]]), axis=1)
                self.status.set('WF 1 and 2 combined')
        else:
            self.status.set('number of WF 0 and WF 1 should be same')
            
    def shift_wfs(self):
        self.shift_wf0 = {}
        self.shift_wf1 = {}
        self.shift_time = {}
        OP_switch = 1
        SO_switch = 1
        try:
            time_delays = self.so_optimal_go.all_time_delays[self.n0.get()]
        except AttributeError:
            OP_switch = 0
        else:
            time_delays = self.so_optimal_go.all_time_delays[self.n0.get()]
            for key, time_delay in zip(self.interpo_wf0.keys(), time_delays):
                newkey = 'OP_shifted_'+key
                # idx0 = np.where(self.t[key]>=self.t_left)[0][0]
                # idx1 = np.where(self.t[key]<=self.t_right)[0][-1]
                self.shift_wf0[key] = self.interpo_wf0[key](self.t[key]-time_delay)
                # self.shift_time[key] = self.t[key]
                self.x[newkey] = self.shift_wf0[key]
                self.t[newkey] = self.t[key]
                zero_filler0 = np.zeros_like(self.x[newkey])
                spec_shift0 = np.stack((self.t[newkey], self.x[newkey], zero_filler0), axis = 1)      
                np.savetxt(self.path + 'OP_shifted_' + key + '_1.dat', spec_shift0)
            for key, time_delay in zip(self.interpo_wf1.keys(), time_delays):
                # idx0 = np.where(self.t[key]>=self.t_left)[0][0]
                # idx1 = np.where(self.t[key]<=self.t_right)[0][-1]
                # self.shift_time[key] = self.t[key]
                newkey = 'OP_shifted_'+ key
                self.shift_wf1[key] = self.interpo_wf1[key](self.t[key]-time_delay)
                self.x[newkey] = self.shift_wf1[key]
                self.t[newkey] = self.t[key]
                zero_filler1 = np.zeros_like(self.x[newkey])
                spec_shift1 = np.stack((self.t[newkey], self.x[newkey], zero_filler1), axis = 1)
                np.savetxt(self.path + 'OP_shifted_' + key + '_1.dat', spec_shift1)
        
        try: 
            time_delays = self.so_phase_go.all_time_delays[self.n0.get()]
        except AttributeError:
            SO_switch = 0
        else:
            time_delays = self.so_phase_go.all_time_delays[self.n0.get()]
            for key, time_delay in zip(self.interpo_wf0.keys(), time_delays):
                newkey = 'SO_shifted_'+ key
                # idx0 = np.where(self.t[key]>=self.t_left)[0][0]
                # idx1 = np.where(self.t[key]<=self.t_right)[0][-1]
                self.shift_wf0[key] = self.interpo_wf0[key](self.t[key]-time_delay)
                # self.shift_time[key] = self.t[key]
                self.x[newkey] = self.shift_wf0[key]
                self.t[newkey] = self.t[key]
                zero_filler0 = np.zeros_like(self.x[newkey])
                spec_shift0 = np.stack((self.t[newkey], self.x[newkey], zero_filler0), axis = 1)      
                np.savetxt(self.path + 'SO_shifted_' + key + '_1.dat', spec_shift0)
            for key, time_delay in zip(self.interpo_wf1.keys(), time_delays):
                # idx0 = np.where(self.t[key]>=self.t_left)[0][0]
                # idx1 = np.where(self.t[key]<=self.t_right)[0][-1]
                newkey = 'SO_shifted_'+ key
                self.shift_wf1[key] = self.interpo_wf1[key](self.t[key]-time_delay)
                # self.shift_time[key] = self.t[key]
                self.x[newkey] = self.shift_wf1[key]
                self.t[newkey] = self.t[key]
                zero_filler1 = np.zeros_like(self.x[newkey])
                spec_shift1 = np.stack((self.t[newkey], self.x[newkey], zero_filler1), axis = 1)
                np.savetxt(self.path + 'SO_shifted_' + key + '_1.dat', spec_shift1)
                
        try: 
            time_delays = self.so_optJ_go.all_time_delays[self.n0.get()]
        except AttributeError:
            OPJ_switch = 0
        else:
            time_delays = self.so_optJ_go.all_time_delays[self.n0.get()]
            for key, time_delay in zip(self.interpo_wf0.keys(), time_delays):
                newkey = 'OPJ_shifted_'+ key
                # idx0 = np.where(self.t[key]>=self.t_left)[0][0]
                # idx1 = np.where(self.t[key]<=self.t_right)[0][-1]
                self.shift_wf0[key] = self.interpo_wf0[key](self.t[key]-time_delay)
                # self.shift_time[key] = self.t[key]
                self.x[newkey] = self.shift_wf0[key]
                self.t[newkey] = self.t[key]
                zero_filler0 = np.zeros_like(self.x[newkey])
                spec_shift0 = np.stack((self.t[newkey], self.x[newkey], zero_filler0), axis = 1)      
                np.savetxt(self.path + 'OPJ_shifted_' + key + '_1.dat', spec_shift0)
            for key, time_delay in zip(self.interpo_wf1.keys(), time_delays):
                # idx0 = np.where(self.t[key]>=self.t_left)[0][0]
                # idx1 = np.where(self.t[key]<=self.t_right)[0][-1]
                newkey = 'OPJ_shifted_'+ key
                self.shift_wf1[key] = self.interpo_wf1[key](self.t[key]-time_delay)
                # self.shift_time[key] = self.t[key]
                self.x[newkey] = self.shift_wf1[key]
                self.t[newkey] = self.t[key]
                zero_filler1 = np.zeros_like(self.x[newkey])
                spec_shift1 = np.stack((self.t[newkey], self.x[newkey], zero_filler1), axis = 1)
                np.savetxt(self.path + 'OPJ_shifted_' + key + '_1.dat', spec_shift1)
                
                
        self.listbox_filenames_content.set(list(self.x))
        self.filenames = list(self.x)
        if SO_switch == 0 and OP_switch == 0 and OPJ_switch == 0:
            self.status.set('Error: No phase information is provided')
        else:
            self.status.set('Shifted pulses calculated')
        self.so_switch = SO_switch
        self.op_switch = OP_switch
        
    # def so_calculator(self):
    #     progress = Queue()
    #     so_simulation = Process(target = self.simulate_so, 
        
    def simulate_so(self):
        try:
            time_delays = self.so_phase_go.all_time_delays
        except AttributeError:
            type_key = 'Interference_'
        else:
            type_key = 'SO_'
            
        newkey = type_key + self.bf.find_list_common(self.filenames_wf_0)
        self.x[newkey] = sum(self.x[key] for key in self.filenames_wf_0)
        self.t[newkey] = self.t[self.filenames_wf_0[0]]
        zero_filler = np.zeros_like(self.x[newkey])
        spec = np.stack((self.t[newkey],self.x[newkey],zero_filler), axis = 1)
        np.savetxt(self.path + newkey + '_1.dat', spec)
        
        newkey = type_key + self.bf.find_list_common(self.filenames_wf_1)
        self.x[newkey] = sum(self.x[key] for key in self.filenames_wf_1)
        self.t[newkey] = self.t[self.filenames_wf_1[0]]
        zero_filler = np.zeros_like(self.x[newkey])
        spec = np.stack((self.t[newkey],self.x[newkey],zero_filler), axis = 1)
        np.savetxt(self.path + newkey + '_1.dat', spec)
        

            
        self.listbox_filenames_content.set(list(self.x))
        self.filenames = list(self.x)
        
        self.status.set('simulated superoscillation calculated')

        

    def start_thread_cal_so(self):
        self.status.set('calculating time delays, please wait...')
        self.window_so_status.update()
        self.progressbar.start()
        newthread = Thread(target = self.calculate_so_phases)
        newthread.start()
        
    def start_thread_maxJ_var_window(self):
        self.J_max = []
        self.localf_ave = []
        self.t_delay_save = []
        self.field = []
        self.status.set('calculating best J for different window, please wait...')
        self.window_so_status.update()
        self.progressbar.start()
        newthread = Thread(target = self.max_J_var_window)
        newthread.start()
        
    def start_thread_cal_optJ(self):
        self.status.set('calculating time delays, please wait...')
        self.window_so_status.update()
        self.progressbar.start()
        newthread = Thread(target = self.calculate_optimal_phase)
        newthread.start()
        
        
    def calculate_so_phases(self):
        # self.status.set('calculating time delays, please wait...')
        # self.window_so_status.update()
        # self.progressbar.start()
        self.so_phase_go = socal.Calculate_SO_phase(self.path, self.filenames_wf_0)
        # self.t_left = self.t[self.filenames_wf_0[0]][0]
        # self.t_right = self.t[self.filenames_wf_0[0]][-1]
        self.progressbar.stop()
        self.status.set('SO phases calculated')
        
    def calculate_optimal_phase(self):
        # self.status.set('calculating time delays, please wait...')
        # self.window_so_status.update()
        self.so_optimal_go = socal.Calculate_optimal_phase(self.path, self.filenames_wf_0, self.filenames_wf_1)
        # self.t_left = self.t[self.filenames_wf_0[0]][0]
        # self.t_right = self.t[self.filenames_wf_0[0]][-1]
        self.progressbar.stop()
        self.status.set('Optimal phases calculated')
        
        
    def max_J_var_window(self):
        half_window_len = np.arange(0.2, 4, step = 0.2)
        for idx, h_win in enumerate(half_window_len):
            self.so_optJ_go = socal.Calculate_optimal_phase_varwindow(self.path, self.filenames_wf_0, self.filenames_wf_1, 
                                                            window = (-h_win, h_win))
            self.progressbar['value'] = float(idx/(len(half_window_len)-1)*100)
            self.J_max.append(self.so_optJ_go.best_J)
            self.localf_ave.append(self.so_optJ_go.localf_ave)
            # self.t_delay_save.append(self.so_optJ_go.all_time_delays[0])
            
            self.window_so_status.update()
        self.J_max = np.array(self.J_max)
        zero_filler = np.zeros_like(self.J_max)
        spec = np.stack((half_window_len, np.log10(self.J_max), zero_filler), axis = 1)
        np.savetxt(self.path+'Jmax_1.dat', spec)
        self.progressbar.stop()
        self.status.set('Best J for different window sizes calculated')
        self.window_so_status.update()
        
    def find_closest_so(self):
        shift_wf0 = {}
        x_all_ref = self.x[self.filenames_so_0[0]]
        t_all_ref = self.t[self.filenames_so_0[0]]
        x_so_ref = self.get_so_region(t_all_ref, x_all_ref)
        max_ref = max(x_so_ref)
        diff0 = sum(abs(x_so_ref-0))
        for i in range(0,500):
            time_delays = self.so_phase_go.all_time_delays[i]
            for key, time_delay in zip(self.interpo_wf0.keys(), time_delays):
                shift_wf0[key] = self.interpo_wf0[key](self.t[key]-time_delay)
                x_all_sam = sum(shift_wf0[key] for key in shift_wf0.keys())
                t_all_sam = self.t[self.filenames_wf_0[0]]
                x_so_sam = self.get_so_region(t_all_sam, x_all_sam)
                max_sam = max(x_so_sam)
                x_so_sam = x_so_sam*(max_ref/max_sam)
                diff1 = sum(abs(x_so_ref - x_so_sam))
                if diff1 < diff0:
                    diff0 = diff1
                    self.n0.set(i)
        
        
    def get_so_region(self, t, x, l_tw = 1.2, t_step = 0.05):
        t_points = round(l_tw/t_step)
        t_start = (t[0] + t[-1])/2 - l_tw/2
        idx_ts = np.where(t>=t_start)[0][0]
        idx_te = idx_ts + t_points
        x_so = x[idx_ts:idx_te]
        return x_so
    
    def discrimin_normal(self):
        timestep = 0.05
        if len(self.idx_wf_0) == len(self.idx_wf_1) and len(self.idx_so_0) == len(self.idx_so_1):
            # calculate discriminability of basic waveforms
            for i in range (0, len(self.idx_wf_0)):
                newkey = 'WF_discriminability_'+ str(self.n_dis)+ '_' + self.bf.find_str_common(self.filenames_wf_0[i], self.filenames_wf_1[i])
                self.J[newkey] = []
                self.J_base[newkey] = []
                self.J_log[newkey] = []
                self.tw_len[newkey] = []
                t = self.t[self.filenames[self.idx_wf_0[i]]]
                E1 = self.x[self.filenames[self.idx_wf_0[i]]]
                E2 = self.x[self.filenames[self.idx_wf_1[i]]]
                self.idx_peaks = find_peaks(abs(E1), height = 0.5*max(E1))[0]
                self.idx_peaks_reduce = self.idx_peaks - round(len(t)/2)
                if len(self.idx_peaks[self.idx_peaks_reduce > 0]) == 0:
                    idx0 = self.idx_peaks[-1]
                else:
                    idx0 = self.idx_peaks[self.idx_peaks_reduce > 0][0]
                idx1 = idx0 + 1
                # if self.timewindow_bot.get() != 0 or self.timewindow_top.get() != 0:
                #     tw_top = self.timewindow_top.get()
                #     tw_bot = self.timewindow_bot.get()
                # else:
                #     tw_bot = t[0]
                #     tw_top = t[-1]
                # idx1 = np.where(t <= tw_top)[0][-1]
                # idx0 = np.where(t >= tw_bot)[0][0]
                while idx0 >= 0 and idx1 <= len(t)-1:
                    t_w = abs(t[idx1]-t[idx0])
                    E1_w = E1[idx0:idx1]
                    E2_w = E2[idx0:idx1]
                    # peak1 = max(abs(E1_w))
                    self.tw_len[newkey].append(t[idx1] - t[idx0])
                    J_temp = 2*simps(abs(E1_w-E2_w)**2)/simps(abs(E1_w)**2+abs(E2_w)**2)
                    # J_temp = simps(abs(E1_w-E2_w)**2)/t_w
                    self.J[newkey].append(J_temp)
                    self.J_base[newkey].append(simps((abs(E1_w-E2_w)**2), dx=timestep)/abs(t[idx1]-t[idx0]))
                    self.J_log[newkey].append(np.log10(J_temp))
                    idx0 = idx0 - 1
                    idx1 = idx1 + 1
                self.tw_len[newkey] = np.array(self.tw_len[newkey])
                self.J[newkey] = np.array(self.J[newkey])
                self.J_base[newkey] = np.array(self.J_base[newkey])
                self.J_log[newkey] = np.array(self.J_log[newkey])
                dp.Basic_functions.save_data(self.path+'/J calculation results/'+newkey+'_1.dat', self.tw_len[newkey], self.J_log[newkey])
            # calculate discriminability of superoscillations
            for j in range(0, len(self.idx_so_0)):
                newkey = 'SO_discriminability_' + str(self.n_dis) + '_' +self.bf.find_str_common(self.filenames_so_0[j], self.filenames_so_1[j])
                self.J[newkey] = []
                self.J_base[newkey] = []
                self.J_log[newkey] = []
                self.tw_len[newkey] = []
                t = self.t[self.filenames[self.idx_so_0[j]]]
                E1 = self.x[self.filenames[self.idx_so_0[j]]]
                E2 = self.x[self.filenames[self.idx_so_1[j]]]
                if self.timewindow_bot.get() != 0 or self.timewindow_top.get() != 0:
                    tw_top = self.timewindow_top.get()
                    tw_bot = self.timewindow_bot.get()
                else:
                    tw_bot = t[0]
                    tw_top = t[-1]
                idx1 = np.where(t<tw_top)[0][-1]
                idx0 = np.where(t>tw_bot)[0][0]
                while idx0 < idx1:
                    t_w = abs(t[idx1]-t[idx0])
                    E1_w = E1[idx0:idx1]
                    E2_w = E2[idx0:idx1]
                    # peak1 = max(abs(E1_w))
                    J_temp = 2*simps(abs(E1_w-E2_w)**2)/simps(abs(E1_w)**2+abs(E2_w)**2)
                    # J_temp = simps(abs(E1_w-E2_w)**2)/t_w
                    self.tw_len[newkey].append(t[idx1] - t[idx0])
                    self.J[newkey].append(J_temp)
                    self.J_base[newkey].append(simps((abs(E1_w-E2_w)**2), dx=timestep)/abs(t[idx1]-t[idx0]))
                    self.J_log[newkey].append(np.log10(J_temp))
                    idx0 = idx0 + 1
                    idx1 = idx1 - 1
                self.tw_len[newkey] = np.array(self.tw_len[newkey])
                self.J[newkey] = np.array(self.J[newkey])
                self.J_base[newkey] = np.array(self.J_base[newkey])
                self.J_log[newkey] = np.array(self.J_log[newkey])
                dp.Basic_functions.save_data(self.path+'/J calculation results/'+newkey+'_1.dat', self.tw_len[newkey], self.J_log[newkey])
        else:
            if len(self.idx_wf_0) != len(self.idx_wf_1):
                self.status.set('Correct Error: select same number of waveform1 and waveform2')
            if len(self.idx_so_0) != len(self.idx_so_1):
                self.status.set('Correct Error: select same number of superosc1 and superosc2')
        # plt.close()        
        # self.canvas0.figure, self.ax = plt.subplots(1, 1, figsize = (16,8))
        # for key in list(self.J):
        #     self.ax.plot(self.tw_len[key], self.J[key], label = key)
        #     self.ax.legend(loc = 'best')
        #     self.ax.set_xlabel('time window size (ps)', fontsize = 15)
        #     self.ax.set_ylabel('discriminability (arb.u.)', fontsize = 15)
        #     self.ax.tick_params(axis='both', labelsize=15)
        # self.canvas0.draw()
        self.n_dis = self.n_dis + 1
        if self.n_dis > 3:
            self.tw_len = {}
            self.J = {}
            self.n_dis = 0
        
        new_window = tk.Toplevel(self.root)
        new_window.title('Discriminability time domain')
        if self.discrimin_type.get() == 'normalized':
            new_window.title('normalized discriminability')
            plot_dis = pt.Plottool(new_window, self.tw_len, self.J)
        if self.discrimin_type.get() == 'original':   
            new_window.title('non normalized discriminability')
            plot_dis = pt.Plottool(new_window, self.tw_len, self.J_base)
        if self.discrimin_type.get() == 'log scale': 
            new_window.title('log scale normalized discriminability')
            plot_dis = pt.Plottool(new_window, self.tw_len, self.J_log)
    
    def discrimin_polar_J(self):
        timestep = 0.05
        if len(self.idx_wf_0) == len(self.idx_wf_1) and len(self.idx_so_0) == len(self.idx_so_1):
            if len(self.idx_wf_0) == len(self.idx_so_0):
            # calculate discriminability of basic waveforms
                for i in range (0, len(self.idx_wf_0)):
                    newkey = 'J_polar_'+ str(self.n_dis)+ '_' + self.bf.find_str_common(self.filenames_wf_0[i], self.filenames_wf_1[i])
                    self.J[newkey] = []
                    self.J_base[newkey] = []
                    self.J_log[newkey] = []
                    self.tw_len[newkey] = []
                    t = self.t[self.filenames[self.idx_wf_0[i]]]
                    E1_x = self.x[self.filenames[self.idx_wf_0[i]]]
                    E1_y = self.x[self.filenames[self.idx_so_0[i]]]
                    E2_x = self.x[self.filenames[self.idx_wf_1[i]]]
                    E2_y = self.x[self.filenames[self.idx_so_1[i]]]
                    # self.idx_peaks = find_peaks(abs(E1), height = 0.5*max(E1))[0]
                    # self.idx_peaks_reduce = self.idx_peaks - round(len(t)/2)
                    # if len(self.idx_peaks[self.idx_peaks_reduce > 0]) == 0:
                    #     idx0 = self.idx_peaks[-1]
                    # else:
                    #     idx0 = self.idx_peaks[self.idx_peaks_reduce > 0][0]
                    
                    idx0 = int(np.floor(len(t)/2))
                    idx1 = idx0 + 1
                    while idx0 >= 0 and idx1 <= len(t)-1:
                        t_w = abs(t[idx1]-t[idx0])
                        E1x_w = E1_x[idx0:idx1]
                        E2x_w = E2_x[idx0:idx1]
                        E1y_w = E1_y[idx0:idx1]
                        E2y_w = E2_y[idx0:idx1]
                        self.tw_len[newkey].append(t[idx1] - t[idx0])
                        J_temp = 2*simps(abs(E1x_w-E2x_w)**2+abs(E1y_w-E2y_w)**2)/simps(abs(E1x_w)**2+abs(E2x_w)**2+E1y_w**2+E2y_w**2)
                        # J_temp = simps(abs(E1_w-E2_w)**2)/t_w
                        self.J[newkey].append(J_temp)
                        # self.J_base[newkey].append(simps((abs(E1_w-E2_w)**2), dx=timestep)/abs(t[idx1]-t[idx0]))
                        self.J_log[newkey].append(np.log10(J_temp))
                        idx0 = idx0 - 1
                        idx1 = idx1 + 1
                    self.tw_len[newkey] = np.array(self.tw_len[newkey])
                    self.J[newkey] = np.array(self.J[newkey])
                    self.J_base[newkey] = np.array(self.J_base[newkey])
                    self.J_log[newkey] = np.array(self.J_log[newkey])
                    dp.Basic_functions.save_data(self.path+'/J calculation results/'+newkey+'_1.dat', self.tw_len[newkey], self.J_log[newkey])
            else:
                self.status.set('Correct Error: number of x and y components needs to match')
        
        else:
            if len(self.idx_wf_0) != len(self.idx_wf_1):
                self.status.set('Correct Error: select same number of x components')
            if len(self.idx_so_0) != len(self.idx_so_1):
                self.status.set('Correct Error: select same number of y components')
        self.n_dis = self.n_dis + 1
        if self.n_dis > 3:
            self.tw_len = {}
            self.J = {}
            self.n_dis = 0
        
        new_window = tk.Toplevel(self.root)
        new_window.title('Discriminability time domain')
        if self.discrimin_type.get() == 'normalized':
            new_window.title('normalized discriminability')
            plot_dis = pt.Plottool(new_window, self.tw_len, self.J)
        if self.discrimin_type.get() == 'original':   
            new_window.title('non normalized discriminability')
            plot_dis = pt.Plottool(new_window, self.tw_len, self.J_base)
        if self.discrimin_type.get() == 'log scale': 
            new_window.title('log scale normalized discriminability')
            plot_dis = pt.Plottool(new_window, self.tw_len, self.J_log)
            
    
        
    def discrimin_nd_shift(self):
        timestep = 0.05
        len_t = round(0.6/0.05)
        if len(self.idx_wf_0) == len(self.idx_wf_1) and len(self.idx_so_0) == len(self.idx_so_1):
            # calculate discriminability of basic waveforms
            for i in range (0, len(self.idx_wf_0)):
                newkey = 'WF_discriminability_'+ str(self.n_dis)+ '_' + self.bf.find_str_common(self.filenames[self.idx_wf_0[i]], self.filenames[self.idx_wf_1[i]])
                self.J[newkey] = []
                self.J_base[newkey] = []
                self.J_log[newkey] = []
                self.tw_len[newkey] = []
                t = self.t[self.filenames[self.idx_wf_0[i]]]
                E1 = self.x[self.filenames[self.idx_wf_0[i]]]
                E2 = self.x[self.filenames[self.idx_wf_1[i]]]
                
                idx0 = np.where(t>t[0]+0.6)[0][0]
                while idx0 + len_t < len(t):
                    t_w = 1.2
                    E1_w = E1[idx0-len_t:idx0 + len_t+1]
                    E2_w = E2[idx0-len_t:idx0 + len_t+1]
                    # peak1 = max(abs(E1_w))
                    J_temp = 2*simps(abs(E1_w-E2_w)**2)/simps(abs(E1_w)**2+abs(E2_w)**2)
                    # J_temp = simps(abs(E1_w-E2_w)**2)/1.2
                    self.tw_len[newkey].append(t[idx0])
                    self.J[newkey].append(J_temp)
                    self.J_base[newkey].append(simps((abs(E1_w**2-E2_w**2)), dx=timestep)/1.2)
                    self.J_log[newkey].append(np.log10(J_temp))
                    idx0 = idx0 + 1
                self.tw_len[newkey] = np.array(self.tw_len[newkey])
                self.J[newkey] = np.array(self.J[newkey])
                self.J_base[newkey] = np.array(self.J_base[newkey])
                self.J_log[newkey] = np.array(self.J_log[newkey])
                dp.Basic_functions.save_data(self.path+'/J calculation results/'+newkey+'_1.dat', self.tw_len[newkey], self.J_log[newkey])
            # calculate discriminability of superoscillations
            for j in range(0, len(self.idx_so_0)):
                newkey = 'Z_SO_discriminability_' + str(self.n_dis) + '_' +self.bf.find_str_common(self.filenames[self.idx_so_0[j]], self.filenames[self.idx_so_1[j]])
                self.J[newkey] = []
                self.J_base[newkey] = []
                self.J_log[newkey] = []
                self.tw_len[newkey] = []
                t = self.t[self.filenames[self.idx_so_0[j]]]
                E1 = self.x[self.filenames[self.idx_so_0[j]]]
                E2 = self.x[self.filenames[self.idx_so_1[j]]]
                idx0 = np.where(t>t[0]+0.6)[0][0]
                while idx0 + len_t < len(t):
                    E1_w = E1[idx0-len_t:idx0 + len_t]
                    E2_w = E2[idx0-len_t:idx0 + len_t]
                    # peak1 = max(abs(E1_w))
                    J_temp = 2*simps(abs(E1_w-E2_w)**2)/simps(abs(E1_w)**2+abs(E2_w)**2)
                    # J_temp = simps(abs(E1_w-E2_w)**2)/1.2
                    self.tw_len[newkey].append(t[idx0])
                    self.J[newkey].append(J_temp)
                    self.J_base[newkey].append(simps((abs(E1_w**2-E2_w**2)), dx=timestep)/1.2)
                    self.J_log[newkey].append(np.log10(J_temp))
                    idx0 = idx0 + 1
                self.tw_len[newkey] = np.array(self.tw_len[newkey])
                self.J[newkey] = np.array(self.J[newkey])
                self.J_base[newkey] = np.array(self.J_base[newkey])
                self.J_log[newkey] = np.array(self.J_log[newkey])
                dp.Basic_functions.save_data(self.path+'/J calculation results/'+newkey+'_1.dat', self.tw_len[newkey], self.J_log[newkey])
        else:
            if len(self.idx_wf_0) != len(self.idx_wf_1):
                self.status.set('Correct Error: select same number of waveform1 and waveform2')
            if len(self.idx_so_0) != len(self.idx_so_1):
                self.status.set('Correct Error: select same number of superosc1 and superosc2')
        
        new_window = tk.Toplevel(self.root)
        new_window.title('Discriminability time domain')
        if self.discrimin_type.get() == 'normalized':
            new_window.title('normalized discriminability')
            plot_dis = pt.Plottool(new_window, self.tw_len, self.J)
        if self.discrimin_type.get() == 'original':   
            new_window.title('non normalized discriminability')
            plot_dis = pt.Plottool(new_window, self.tw_len, self.J_base)
        if self.discrimin_type.get() == 'log scale': 
            new_window.title('log scale normalized discriminability')
            plot_dis = pt.Plottool(new_window, self.tw_len, self.J_log)
            
        self.n_dis = self.n_dis + 1
        if self.n_dis > 3:
            self.tw_len = {}
            self.J = {}
            self.n_dis = 0
            
    def discrimin_normal_fd(self):
        # fstep = 0.
        # dE_sum1 = 0
        # dE_sum2 = 0
        self.f, self.sx = dp.Basic_functions().fftx(self.x, self.t, 2)
        self.Ef = dp.Basic_functions().dict_getabs(self.sx)
        if len(self.idx_wf_0) == len(self.idx_wf_1) and len(self.idx_so_0) == len(self.idx_so_1):
            # calculate discriminability of basic waveforms
            for i in range (0, len(self.idx_wf_0)):
                newkey = 'WF_discriminability_'+ str(self.n_dis)+ '_' + self.bf.find_str_common(self.filenames[self.idx_wf_0[i]], self.filenames[self.idx_wf_1[i]])
                
                f = self.f[self.filenames[self.idx_wf_0[i]]]
                # sx1 = self.sx[self.filenames[self.idx_wf_0[i]]]
                # sx2 = self.sx[self.filenames[self.idx_wf_1[i]]]
                E1 = abs(self.sx[self.filenames[self.idx_wf_0[i]]])
                E2 = abs(self.sx[self.filenames[self.idx_wf_1[i]]])
                idx1 = np.where(f<1.6)[0][-1]
                idx0 = np.where(f>0.1)[0][0]
                # sx1_w = sx1[idx0:idx1]
                # sx2_w = sx2[idx0:idx1]
                E1_w = E1[idx0:idx1]
                E2_w = E2[idx0:idx1]
                # dE_sum1 = dE_sum1 + sx1_w
                # dE_sum2 = dE_sum2 + sx2_w
                # peak1 = max([max(E1_w),max(E2_w)])
                self.f_J[newkey] = f[idx0:idx1]
                # self.J[newkey] = abs(E1_w-E2_w)/peak1
                # self.J_base[newkey] = abs(E1_w-E2_w)/peak1
                # self.J_log[newkey] = np.log10(abs(E1_w-E2_w)/peak1)  
                self.J[newkey] = abs(E1_w-E2_w)
                self.J_base[newkey] = abs(E1_w-E2_w)
                self.J_log[newkey] = np.log10(abs(E1_w-E2_w))   
            # calculate discriminability of superoscillations
            # newkey2 = 'WF_discriminability_sum_' + str(self.n_dis)
            # self.J[newkey2] = abs(dE_sum1 - dE_sum2)
            # self.f_J[newkey2] = f[idx0:idx1]
            for j in range(0, len(self.idx_so_0)):
                newkey = 'SO_discriminability_' + str(self.n_dis) + '_' +self.bf.find_str_common(self.filenames[self.idx_so_0[j]], self.filenames[self.idx_so_1[j]])
                f = self.f[self.filenames[self.idx_so_0[j]]]
                E1 = abs(self.sx[self.filenames[self.idx_so_0[j]]])
                E2 = abs(self.sx[self.filenames[self.idx_so_1[j]]])
                idx1 = np.where(f<1.6)[0][-1]
                idx0 = np.where(f>0.1)[0][0]
                E1_w = E1[idx0:idx1]
                E2_w = E2[idx0:idx1]
                # peak1 = max([max(E1_w),max(E2_w)])
                self.f_J[newkey] = f[idx0:idx1]
                # self.J[newkey] = abs(E1_w-E2_w)/peak1
                # self.J_base[newkey] = abs(E1_w-E2_w)/peak1
                # self.J_log[newkey] = np.log10(abs(E1_w-E2_w)/peak1)
                self.J[newkey] = abs(E1_w-E2_w)
                self.J_base[newkey] = abs(E1_w-E2_w)
                self.J_log[newkey] = np.log10(abs(E1_w-E2_w))
        else:
       
            if len(self.idx_wf_0) != len(self.idx_wf_1):
                self.status.set('Correct Error: select same number of waveform1 and waveform2')
                if len(self.idx_so_0) != len(self.idx_so_1):
                    self.status.set('Correct Error: select same number of waveform1 and waveform2' + '\n'
                                    +'Correct Error: select same number of superosc1 and superosc2')
            else:
                if len(self.idx_so_0) != len(self.idx_so_1):
                    self.status.set('Correct Error: select same number of superosc1 and superosc2')

        
        new_window = tk.Toplevel(self.root)
        new_window2 = tk.Toplevel(self.root)
        new_window.title('Frequency Domain Discriminability')
        if self.discrimin_type.get() == 'normalized':
            new_window.title('normalized discriminability')
            plot_fd = pt.Plottool(new_window2, self.f, self.Ef)
            plot_dis = pt.Plottool(new_window, self.f_J, self.J)
        if self.discrimin_type.get() == 'original':   
            new_window.title('non normalized discriminability')
            plot_dis = pt.Plottool(new_window, self.f_J, self.J_base)
        if self.discrimin_type.get() == 'log scale': 
            new_window.title('log scale normalized discriminability')
            plot_dis = pt.Plottool(new_window, self.f_J, self.J_log)        
        self.n_dis = self.n_dis + 1
        if self.n_dis > 3:
            self.tw_len = {}
            self.J = {}
            self.J_base = {}
            self.J_log = {}
            self.f_J = {}
            self.n_dis = 0
        
   
    
    
        
    
    def D_discrimin(self):
        D_all = []
        for i in range(0,self.n_dis):
            D_wf = 0
            D_so = 0
            n_wf = 0
            n_so = 0
            self.i = i
            for key in self.tw_len.keys():
                if 'SO_discriminability_'+str(i) in key:
                    t = self.tw_len[key][np.where(self.tw_len[key]<=1.2)]
                    J = self.J[key][np.where(self.tw_len[key]<=1.2)]
                    D_so = D_so + simps(J, dx = abs(t[0]-t[1]))
                    n_so = n_so + 1
                    self.n_so = n_so
                    
                if 'WF_discriminability_'+str(i) in key:
                    t = self.tw_len[key][np.where(self.tw_len[key]<=1.2)]
                    J = self.J[key][np.where(self.tw_len[key]<=1.2)]
                    D_wf = simps(J, dx = abs(t[0]-t[1])) + D_wf
                    n_wf = n_wf + 1
                    self.n_wf = n_wf
            D_all.append(np.log10(D_so/n_so/(D_wf/n_wf)))
                
        self.D_J.set(str(D_all))
            
    def fftx_stft(self, t, x, pad = 2):
        ts = abs(t[1]-t[0])
        NFFT = int(2**(np.ceil(np.log2(len(t)))+pad))
        fs = 1/ts
        x_normal = x-np.average(x)
        # freq_auto = 0
        # t_auto = 0
        # zx = 0
        freq_auto, t_auto, zx = stft(x_normal, fs, nperseg = len(t)/3, noverlap=len(t)/3-1,
                                            window = ('gaussian', 11), nfft = NFFT)
        return freq_auto, t_auto, zx  
    
    def Localf_stft(self):
        newwindow = tk.Toplevel()
        
        if len(self.idx_wf_0) != 0:
            for i in self.idx_wf_0:
                key = self.filenames[i]
                newkey = 'localf_'+key
                # E_t = self.x[key]
                freq, time, zx = self.fftx_stft(self.t[key], self.x[key])
                f_ave = np.zeros_like(time)
                nrow, ncol = np.shape(zx)
                ampl = abs(zx)
                for i in range(0, ncol):
                    ampl_sum = sum(ampl[:,i])
                    f_ave[i] = sum(freq*ampl[:,i]/ampl_sum)
                self.t[newkey] = self.t[key][0:len(time)-1]
                self.x[newkey] = f_ave[0:len(time)-1]
                
        if len(self.idx_wf_1) != 0:
            for i in self.idx_wf_1:
                key = self.filenames[i]
                newkey = 'localf_'+key
                # E_t = self.x[key]
                freq, time, zx = self.fftx_stft(self.t[key], self.x[key])
                f_ave = np.zeros_like(time)
                nrow, ncol = np.shape(zx)
                ampl = abs(zx)
                for i in range(0, ncol):
                    ampl_sum = sum(ampl[:,i])
                    f_ave[i] = sum(freq*ampl[:,i]/ampl_sum)
                self.t[newkey] = self.t[key][0:len(time)-1]
                self.x[newkey] = f_ave[0:len(time)-1]
        
        if len(self.idx_so_0) != 0:
            for i in self.idx_so_0:
                key = self.filenames[i]
                newkey = 'localf_'+key
                # E_t = self.x[key]
                freq, time, zx = self.fftx_stft(self.t[key], self.x[key])
                f_ave = np.zeros_like(time)
                nrow, ncol = np.shape(zx)
                ampl = abs(zx)
                for i in range(0, ncol):
                    ampl_sum = sum(ampl[:,i])
                    f_ave[i] = sum(freq*ampl[:,i]/ampl_sum)
                self.t[newkey] = self.t[key][0:len(time)-1]
                self.x[newkey] = f_ave[0:len(time)-1]
        if len(self.idx_so_1) != 0:
            for j in self.idx_so_1:
               key = self.filenames[j]
               newkey = 'localf_'+key
               freq, time, zx = self.fftx_stft(self.t[key], self.x[key])
               f_ave = np.zeros_like(time)
               nrow, ncol = np.shape(zx)
               ampl = abs(zx)
               for i in range(0, ncol):
                   ampl_sum = sum(ampl[:,i])
                   f_ave[i] = sum(freq*ampl[:,i]/ampl_sum)
               self.t[newkey] = self.t[key][0:len(time)-1]
               self.x[newkey] = f_ave[0:len(time)-1]
        plot_localf = pt.Plottool(newwindow, self.t, self.x)
    
                
    def Local_f(self):
        newwindow = tk.Toplevel()
        if len(self.idx_so_0) != 0:
            for i in self.idx_so_0:
                key = self.filenames[i]
                newkey = 'localf_'+key
                # E_t = self.x[key]
                E_cx = hilbert(self.x[key])
                t = self.t[key]
                phase = np.unwrap(np.angle(E_cx))+10*pi
                self.x[newkey] = np.gradient(phase, t[1]-t[0])/2/pi
                # self.x[newkey] = self.locfreq2(t, E_t)
                self.t[newkey] = t
        if len(self.idx_so_1) != 0:
            for j in self.idx_so_1:
                key = self.filenames[j]
                newkey = 'localf_'+key
                E_cx = hilbert(self.x[key])
                # E_t = self.x[key]
                t = self.t[key]
                phase = np.unwrap(np.angle(E_cx))+10*pi
                self.x[newkey] = np.gradient(phase, t[1]-t[0])/2/pi
                # self.x[newkey] = self.locfreq2(t, E_t)
                self.t[newkey] = t
        plot_localf = pt.Plottool(newwindow, self.t, self.x)
        
    def locfreq2(self, t, E):
        
        return np.sqrt(np.abs(
            np.gradient(np.gradient(E, t[1] - t[0]), t[1] - t[0]) / E
        ))
    
    def display_timedelays(self):
        self.td1.set(str(self.so_phase_go.all_time_delays[self.n0.get()]))
        self.td2.set(str(self.so_phase_go.all_time_delays[self.n0.get()+1]))
        self.td3.set(str(self.so_phase_go.all_time_delays[self.n0.get()+2]))
         
    def calculate_phase(self):
        self.ph0 = np.array([self.phase_05.get(),self.phase_06.get(),self.phase_07.get(),self.phase_08.get()])
        try :
            self.ph1.set(str(np.array(self.so_phase_go.all_time_delays[self.n0.get()])+self.ph0))
            self.ph2.set(str(np.array(self.so_phase_go.all_time_delays[self.n0.get()+1])+self.ph0))
            self.ph3.set(str(np.array(self.so_phase_go.all_time_delays[self.n0.get()+2])+self.ph0)) 
        except:
            try:
                self.ph1.set(str(np.array(self.so_optimal_go.all_time_delays[self.n0.get()])+self.ph0))
                self.ph2.set(str(np.array(self.so_optimal_go.all_time_delays[self.n0.get()+1])+self.ph0))
                self.ph3.set(str(np.array(self.so_optimal_go.all_time_delays[self.n0.get()+2])+self.ph0))
            except:
                self.label_status.set('no result to show')
            else:
                self.ph1.set(str(np.array(self.so_optimal_go.all_time_delays[self.n0.get()])+self.ph0))
                self.ph2.set(str(np.array(self.so_optimal_go.all_time_delays[self.n0.get()+1])+self.ph0))
                self.ph3.set(str(np.array(self.so_optimal_go.all_time_delays[self.n0.get()+2])+self.ph0))
        else:
            self.ph1.set(str(np.array(self.so_phase_go.all_time_delays[self.n0.get()])+self.ph0))
            self.ph2.set(str(np.array(self.so_phase_go.all_time_delays[self.n0.get()+1])+self.ph0))
            self.ph3.set(str(np.array(self.so_phase_go.all_time_delays[self.n0.get()+2])+self.ph0)) 
    def plot_all(self):
        newwindow = tk.Toplevel()
        self.so_plotgo = pt.Plottool(newwindow, self.t, self.x)
        
    def plot_predictions(self):
        try:
            all_time_delays = self.so_phase_go.all_time_delays[self.n0.get():self.n0.get()+3]
        except AttributeError:
            try:
                all_time_delays = self.so_optimal_go.all_time_delays[self.n0.get():self.n0.get()+3]
            except AttributeError:
                self.status.set('No result to show')
            else:
                self.plot_predictions_OP()
        else:
            self.plot_predictions_SO()
    
    def plot_predictions_SO(self):
        all_time_delays = self.so_phase_go.all_time_delays[self.n0.get():self.n0.get()+3]
        plot0 = self.ax[0]
        plot1 = self.ax[1]
        plot0.cla()
        plot1.cla()
        plot0.set_title('Best superoscilations')
    
        #long_time_window = 10 * time_window
        long_time_window = self.so_phase_go.pulses[self.filenames[self.idx_wf_0[0]]].time
        for n, time_delays in enumerate(all_time_delays[:3]):
            #plt.plot(long_time_window, get_combined_field(time_delays, long_time_window), label="theory {}".format(n))
            E_c = self.so_phase_go.get_combined_field(time_delays, long_time_window)
            plot0.errorbar(
                long_time_window, 
                E_c, 
                yerr = self.so_phase_go.get_err_combined_field(time_delays, long_time_window),
                label="theory {}".format(n)
            )
            zero_filler = np.zeros_like(long_time_window)
            spec_save = np.stack((long_time_window, E_c, zero_filler), axis = 1)
            np.savetxt(self.path+'best_'+str(n)+'_so_1.dat', spec_save)
       
        
    
        plot0.set_xlabel('time (ps)')
        plot0.set_ylabel('filed (arb. units)')
        plot0.legend()

    
        for n, time_delays in enumerate(all_time_delays[:3]):
            E_c_short = self.so_phase_go.get_combined_field(time_delays, self.so_phase_go.time_window_raw)
            # line = plot1.plot(
            #     self.so_phase_go.time_window_r, self.so_phase_go.get_combined_field(time_delays, self.so_phase_go.time_window), label="theory {}".format(n))
            # line = line[0]
    
            # display raw points
            
            plot1.errorbar(
                self.so_phase_go.time_window_raw, 
                E_c_short, 
                yerr = self.so_phase_go.get_err_combined_field(time_delays, self.so_phase_go.time_window_raw),
                # color = line.get_color(),
                marker='*',
                linestyle='solid',
                label="theory {}".format(n)
            )
            
            # zero_filler = np.zeros_like(self.so_phase_go.time_window_raw)
            # spec_save = np.stack((self.so_phase_go.time_window_raw, E_c_short, zero_filler),
            #                      axis = 1)
            # np.savetxt(self.path+'best_'+str(n)+'_so_short_1.dat', spec_save)
            dp.Basic_functions.save_data(self.path+'best_'+str(n)+'_so_short_1.dat', self.so_phase_go.time_window_raw, E_c_short)
        
        E_hf = self.so_phase_go.pulses[self.so_phase_go.largest_freq].interp_field(self.so_phase_go.time_window_raw)
        plot1.twinx().plot(
            self.so_phase_go.time_window_raw, 
            E_hf, 
            # 'k-',
            color = 'black',
            label=self.so_phase_go.largest_freq,
        )
        
        dp.Basic_functions.save_data(self.path+self.so_phase_go.largest_freq+'_short_1.dat', self.so_phase_go.time_window_raw, E_hf)
       
        plot1.set_xlabel('time (ps)')
        plot1.set_ylabel('filed (arb. units)')
        plot1.legend()
        # plot1.twinx().legend()
    
        self.canvas0.figure.set_tight_layout(True)
        self.canvas0.draw()
    
    def plot_predictions_OP(self):
        all_time_delays = self.so_optimal_go.all_time_delays[self.n0.get():self.n0.get()+3]
        plot0 = self.ax[0]
        plot1 = self.ax[1]
        plot0.cla()
        plot1.cla()
        plot0.set_title('Best superoscilations')
    
        #long_time_window = 10 * time_window
        long_time_window = self.so_optimal_go.pulses0[self.filenames[self.idx_wf_0[0]]].time
        
       
    
        for n, time_delays in enumerate(all_time_delays[:3]):
            #plt.plot(long_time_window, get_combined_field(time_delays, long_time_window), label="theory {}".format(n))
            E1, E2 = self.so_optimal_go.get_combined_fields(time_delays, long_time_window)
            plot0.errorbar(
                long_time_window, 
                E1, 
                yerr = self.so_optimal_go.get_err_combined_field(time_delays, long_time_window),
                label="field 1 theory {}".format(n)
            )
            plot0.errorbar(
                long_time_window, 
                E2, 
                yerr = self.so_optimal_go.get_err_combined_field(time_delays, long_time_window),
                label="field 2 theory {}".format(n)
            )
            zero_filler = np.zeros_like(long_time_window)
            spec_save = np.stack((long_time_window, E1, zero_filler), axis = 1)
            np.savetxt(self.path+'E1_best_'+str(n)+'_op_1.dat', spec_save)
            spec_save = np.stack((long_time_window, E2, zero_filler), axis = 1)
            np.savetxt(self.path+'E2_best_'+str(n)+'_op_1.dat', spec_save)
    
        plot0.set_xlabel('time (ps)')
        plot0.set_ylabel('filed (arb. units)')
        plot0.legend()

    
        for n, time_delays in enumerate(all_time_delays[:3]):
            E1, E2 = self.so_optimal_go.get_combined_fields(time_delays, self.so_optimal_go.time_window)
            line = plot1.plot(
                self.so_optimal_go.time_window, E1, label="field 1 theory {}".format(n))
            line = plot1.plot(
                self.so_optimal_go.time_window, E2, label="field 2 theory {}".format(n))
            line = line[0]
    
            # display raw points
            E1_raw, E2_raw = self.so_optimal_go.get_combined_fields(time_delays, self.so_optimal_go.time_window_raw)
            plot1.errorbar(
                self.so_optimal_go.time_window_raw, 
                E1_raw, 
                yerr = self.so_optimal_go.get_err_combined_field(time_delays, self.so_optimal_go.time_window_raw),
                color = line.get_color(),
                marker='*',
                linestyle=' ',
                
            )
            plot1.errorbar(
                self.so_optimal_go.time_window_raw, 
                E2_raw, 
                yerr = self.so_optimal_go.get_err_combined_field(time_delays, self.so_optimal_go.time_window_raw),
                color = line.get_color(),
                marker='*',
                linestyle=' ',
            )
            
    
        plot1.twinx().plot(
            self.so_optimal_go.time_window, 
            self.so_optimal_go.pulses1[self.so_optimal_go.largest_freq_1].interp_field(self.so_optimal_go.time_window), 
            'k-',
            label=self.so_optimal_go.largest_freq_1,
        )
       
        plot1.set_xlabel('time (ps)')
        plot1.set_ylabel('filed (arb. units)')
        plot1.legend()
    
        self.canvas0.figure.set_tight_layout(True)
        self.canvas0.draw()
        
    def J_all_optimal(self):
        new_window = tk.Toplevel()
        new_window.title('All discriminability')
        window1 = ttk.Labelframe(new_window, text='plot')
        window1.grid(column=0, row=0, sticky='w')
        figure, ax = plt.subplots(1, 1, figsize = (10, 5), dpi = 100)
        figure.set_tight_layout(True)
        canvas0 = FigureCanvasTkAgg(figure = figure, master= window1)
        canvas0.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        long_time_window = self.so_optimal_go.time
        E1, E2 = self.so_optimal_go.get_combined_fields(self.so_optimal_go.all_time_delays[0], long_time_window)
        tw_len, J = dp.Basic_functions.array_cal_J(long_time_window, E1, E2)
        ax.plot(tw_len, np.log10(J), color = 'black', linewidth = 2)
        
        for time_delays in self.so_optimal_go.all_time_delays[self.n0.get():self.n0.get()+100]:
            E1, E2 = self.so_optimal_go.get_combined_fields(time_delays, long_time_window)
            tw_len, J = dp.Basic_functions.array_cal_J(long_time_window, E1, E2)
            ax.plot(tw_len, np.log10(J), color = 'cyan', alpha = 0.05)
            ax.set_xlim((0.2,8))
            ax.set_xlabel('2$T_{obs}$', fontsize = 20)
            ax.set_ylabel('$\log(~J)$', fontsize = 20)
            ax.tick_params(axis='both', labelsize = 20)
        canvas0.draw()   
        toolbar = NavigationToolbar2Tk(canvas0, window1)
            

    def compare_so_op(self):
        try :
            time_delays_0 = self.so_optimal_go.all_time_delays[self.n0.get()]
            time_delays_1 = self.so_phase_go.all_time_delays[self.n0.get()]
        except AttributeError:
            self.status.set('Either SO phases or OP phases are not calculated')
        else:
            plot0 = self.ax[0]
            plot1 = self.ax[1]
            plot0.cla()
            plot1.cla()
            plot0.set_title('Best superoscilations')
        
            #long_time_window = 10 * time_window
            long_time_window = self.so_phase_go.pulses[self.filenames[self.idx_wf_0[0]]].time
            
            short_time_window = self.so_phase_go.time_window
            
        
            E1, E2 = self.so_optimal_go.get_combined_fields(time_delays_0, long_time_window)
            E3 = self.so_phase_go.get_combined_field(time_delays_1, long_time_window)
            E1_s, E2_s = self.so_optimal_go.get_combined_fields(time_delays_0, short_time_window)
            E3_s = self.so_phase_go.get_combined_field(time_delays_1, short_time_window)
            #plt.plot(long_time_window, get_combined_field(time_delays, long_time_window), label="theory {}".format(n))
            
    
            
            plot0.errorbar(
                long_time_window, 
                E1, 
                yerr = self.so_optimal_go.get_err_combined_field(time_delays_0, long_time_window),
                label="OP1 theory"
            )
            plot0.errorbar(
                long_time_window, 
                E2, 
                yerr = self.so_phase_go.get_err_combined_field(time_delays_0, long_time_window),
                label="OP2 theory {}"
            )
            plot0.errorbar(
                long_time_window, 
                E3, 
                yerr = self.so_phase_go.get_err_combined_field(time_delays_1, long_time_window),
                label="SO theory {}"
            )
                
    
        
            plot0.set_xlabel('time (ps)')
            plot0.set_ylabel('filed (arb. units)')
            plot0.legend()
    
        
       
            plot1.plot(
                short_time_window, 
                E1_s, 
                label="OP1 theory"
            )
            plot1.plot(
                short_time_window, 
                E2_s, 
                label="OP2 theory {}"
            )
            plot1.plot(
                short_time_window, 
                E3_s, 
                label="SO theory {}"
            )
           
            plot1.set_xlabel('time (ps)')
            plot1.set_ylabel('filed (arb. units)')
            plot1.legend()
        
            self.canvas0.figure.set_tight_layout(True)
            self.canvas0.draw()
    

            
            
            
