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
import plottool_v5 as pt
import copy

from scipy import stats
from scipy.signal import find_peaks, hilbert, windows
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
        self.filenames = list(t_avg) # get the dict keys of x and t 
        self.status = tk.StringVar(value = 'ready') #s tatus bar textvaribale for displaying status of the application
        self.listbox_filenames_content = tk.StringVar(value = self.filenames) 
        self.J = {}
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
            open('so_tools.ini')
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
        button_calculate_so = ttk.Button(self.window_so_calculate, text = 'Calculate SO', command = self.calculate_so_phases)
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
        button_calculate_so.grid(column = 1, row = 6)
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
        label_D_discrim.grid(column = 0, row = 2, sticky = 'w')
        label_D_discrim_value.grid(column = 1, row = 2, sticky = 'w')
        
        
        '''window content for window_so_plot'''
        
        fig0 = Figure(figsize=(16,8))
        self.canvas0 = FigureCanvasTkAgg(figure=fig0,master=self.window_so_plot)
        self.canvas0.get_tk_widget().grid(column=0,row=0,columnspan=6,rowspan=6,sticky='nsew')
        self.canvas0.figure, self.ax = plt.subplots(2,1,figsize=(16,8))
        self.canvas0.draw()
        
        '''window content for window_so_status'''
        self.label_status = ttk.Label(self.window_so_status, textvariable = self.status, font=('Times', 15))
        
        
        self.label_status.grid(column = 0, row = 0, columnspan = 3)
    
        
        
    
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
        time_delays = self.so_phase_go.all_time_delays[self.n0.get()]
        for key, time_delay in zip(self.interpo_wf0.keys(), time_delays):
            newkey = 'shifted_'+key
            # idx0 = np.where(self.t[key]>=self.t_left)[0][0]
            # idx1 = np.where(self.t[key]<=self.t_right)[0][-1]
            self.shift_wf0[key] = self.interpo_wf0[key](self.t[key]-time_delay)
            # self.shift_time[key] = self.t[key]
            self.x[newkey] = self.shift_wf0[key]
            self.t[newkey] = self.t[key]
            zero_filler0 = np.zeros_like(self.x[newkey])
            spec_shift0 = np.stack((self.t[newkey], self.x[newkey], zero_filler0), axis = 1)      
            np.savetxt(self.path + 'shifted_' + key + '_1.dat', spec_shift0)
        for key, time_delay in zip(self.interpo_wf1.keys(), time_delays):
            # idx0 = np.where(self.t[key]>=self.t_left)[0][0]
            # idx1 = np.where(self.t[key]<=self.t_right)[0][-1]
            self.shift_wf1[key] = self.interpo_wf1[key](self.t[key]-time_delay)
            # self.shift_time[key] = self.t[key]
            newkey = 'shifted_'+ key
            self.x[newkey] = self.shift_wf1[key]
            self.t[newkey] = self.t[key]
            zero_filler1 = np.zeros_like(self.x[newkey])
            spec_shift1 = np.stack((self.t[newkey], self.x[newkey], zero_filler1), axis = 1)
            np.savetxt(self.path + 'shifted_' + key + '_1.dat', spec_shift1)
            
        self.status.set('waveform 1 and 2 shifted')
        
    # def sum_wfs(self):
    #     self.x['sum_wf0'] = sum(self.x[key] for key in self.filenames_wf_0)
    #     self.x['sum_wf1'] = sum(self.x[key] for key in self.filenames_wf_1)
    #     self.t['sum_wf0']
    #     for key1, key2 in zip(self.filenames_wf_0, self.filenames_wf_1):
    #         wf0_sum = 
    def simulate_so(self):
        # self.simu_so0 = {}
        # self.simu_so1 = {}
        newkey_so0 = 'SO_simulated_' + self.bf.find_str_common(self.filenames[self.idx_wf_0[0]], self.filenames[self.idx_wf_0[1]])
        self.x[newkey_so0] = sum(self.shift_wf0[key] for key in self.shift_wf0.keys())
        self.t[newkey_so0] = self.t[self.filenames_wf_0[0]]
        newkey_so1 = 'SO_simulated_' + self.bf.find_str_common(self.filenames[self.idx_wf_1[0]], self.filenames[self.idx_wf_1[1]])
        self.x[newkey_so1] = sum(self.shift_wf1[key] for key in self.shift_wf1.keys())   
        self.t[newkey_so1] = self.t[self.filenames_wf_1[0]]
        # self.selection_so_0.set(newkey_so0)
        # self.selection_so_1.set(newkey_so1)
        self.listbox_filenames_content.set(list(self.x))
        self.filenames = list(self.x)
        zero_filler_1 = np.zeros_like(self.x[newkey_so0])
        zero_filler_2 = np.zeros_like(self.x[newkey_so1])
        spec1 = np.stack((self.t[newkey_so0],self.x[newkey_so0],zero_filler_1), axis = 1)
        spec2 = np.stack((self.t[newkey_so1],self.x[newkey_so1],zero_filler_2), axis = 1)
        np.savetxt(self.path + newkey_so0 + '_1.dat', spec1)
        np.savetxt(self.path + newkey_so1 + '_1.dat', spec2)
        
        self.status.set('simulated superoscillation calculated')
        
        
    def calculate_so_phases(self):
        self.status.set('calculating time delays, please wait...')
        self.window_so_status.update()
        self.so_phase_go = Calculate_SO_phase(self.path, self.filenames_wf_0)
        self.t_left = self.t[self.filenames_wf_0[0]][0]
        self.t_right = self.t[self.filenames_wf_0[0]][-1]
        self.status.set('SO phases calculated')
        
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
                newkey = 'WF_discriminability_'+ str(self.n_dis)+ '_' + self.bf.find_str_common(self.filenames[self.idx_wf_0[i]], self.filenames[self.idx_wf_1[i]])
                self.J[newkey] = []
                self.J_base[newkey] = []
                self.J_log[newkey] = []
                self.tw_len[newkey] = []
                t = self.t[self.filenames[self.idx_wf_0[i]]]
                E1 = self.x[self.filenames[self.idx_wf_0[i]]]
                E2 = self.x[self.filenames[self.idx_wf_1[i]]]
                if self.timewindow_bot.get() != 0 or self.timewindow_top.get() != 0:
                    tw_top = self.timewindow_top.get()
                    tw_bot = self.timewindow_bot.get()
                else:
                    tw_bot = t[0]
                    tw_top = t[-1]
                idx1 = np.where(t<tw_top)[0][-1]
                idx0 = np.where(t>tw_bot)[0][0]
                while idx0 < idx1:
                    E1_w = E1[idx0:idx1]
                    E2_w = E2[idx0:idx1]
                    # peak1 = max(abs(E1_w))
                    self.tw_len[newkey].append(t[idx1] - t[idx0])
                    J_temp = 2*simps(abs(E1_w-E2_w)**2)/simps(abs(E1_w)**2+abs(E2_w)**2)
                    self.J[newkey].append(J_temp)
                    self.J_base[newkey].append(simps((abs(E1_w**2-E2_w**2)), dx=timestep)/abs(t[idx1]-t[idx0]))
                    self.J_log[newkey].append(np.log10(J_temp))
                    idx0 = idx0 + 1
                    idx1 = idx1 - 1
                self.tw_len[newkey] = np.array(self.tw_len[newkey])
                self.J[newkey] = np.array(self.J[newkey])
                self.J_base[newkey] = np.array(self.J_base[newkey])
                self.J_log[newkey] = np.array(self.J_log[newkey])
            # calculate discriminability of superoscillations
            for j in range(0, len(self.idx_so_0)):
                newkey = 'SO_discriminability_' + str(self.n_dis) + '_' +self.bf.find_str_common(self.filenames[self.idx_so_0[j]], self.filenames[self.idx_so_1[j]])
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
                    E1_w = E1[idx0:idx1]
                    E2_w = E2[idx0:idx1]
                    # peak1 = max(abs(E1_w))
                    J_temp = 2*simps(abs(E1_w-E2_w)**2)/simps(abs(E1_w)**2+abs(E2_w)**2)
                    self.tw_len[newkey].append(t[idx1] - t[idx0])
                    self.J[newkey].append(J_temp)
                    self.J_base[newkey].append(simps((abs(E1_w**2-E2_w**2)), dx=timestep)/abs(t[idx1]-t[idx0]))
                    self.J_log[newkey].append(np.log10(J_temp))
                    idx0 = idx0 + 1
                    idx1 = idx1 - 1
                self.tw_len[newkey] = np.array(self.tw_len[newkey])
                self.J[newkey] = np.array(self.J[newkey])
                self.J_base[newkey] = np.array(self.J_base[newkey])
                self.J_log[newkey] = np.array(self.J_log[newkey])
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
            
            
    # def discrimin_convolution(self):
    #     t_ob = np.zeros(int(1.2/0.05)) + 1
    #     if len(self.idx_wf_0) == len(self.idx_wf_1) and len(self.idx_so_0) == len(self.idx_so_1):
    #         # calculate discriminability of basic waveforms
    #         for i in range (0, len(self.idx_wf_0)):
    #             newkey = 'WF_discriminability_'+ str(self.n_dis)+ '_' + self.bf.find_str_common(self.filenames[self.idx_wf_0[i]], 
    #                                                                                             self.filenames[self.idx_wf_1[i]])
    #             t = self.t[self.filenames[self.idx_wf_0[i]]]
    #             E1 = self.x[self.filenames[self.idx_wf_0[i]]]
    #             E2 = self.x[self.filenames[self.idx_wf_1[i]]]
    #             I_diff = abs(E1**2-E2**2)
    #             self.J_conv[newkey] = np.convolve(t_ob, I_diff, mode = 'valid')
    #             self.tw_len[newkey] = t[0:len(self.J_conv[newkey])]
    #         # calculate discriminability of superoscillations
    #         for j in range(0, len(self.idx_so_0)):
    #             newkey = 'SO_discriminability_' + str(self.n_dis) + '_' +self.bf.find_str_common(self.filenames[self.idx_so_0[j]], 
    #                                                                                              self.filenames[self.idx_so_1[j]])
    #             t = self.t[self.filenames[self.idx_so_0[j]]]
    #             E1 = self.x[self.filenames[self.idx_so_0[j]]]
    #             E2 = self.x[self.filenames[self.idx_so_1[j]]]
    #             I_diff = abs(E1**2-E2**2)
    #             self.J_conv[newkey] = np.convolve(t_ob, I_diff, mode = 'valid')
    #             self.tw_len[newkey] = t[0:len(self.J_conv[newkey])]
                
    #     else:
    #         if len(self.idx_wf_0) != len(self.idx_wf_1):
    #             self.status.set('Correct Error: select same number of waveform1 and waveform2')
    #         if len(self.idx_so_0) != len(self.idx_so_1):
    #             self.status.set('Correct Error: select same number of superosc1 and superosc2')
    #     # plt.close()        
    #     # self.canvas0.figure, self.ax = plt.subplots(1, 1, figsize = (16,8))
    #     # for key in list(self.J):
    #     #     self.ax.plot(self.tw_len[key], self.J[key], label = key)
    #     #     self.ax.legend(loc = 'best')
    #     #     self.ax.set_xlabel('time window size (ps)', fontsize = 15)
    #     #     self.ax.set_ylabel('discriminability (arb.u.)', fontsize = 15)
    #     #     self.ax.tick_params(axis='both', labelsize=15)
    #     # self.canvas0.draw()
    #     self.n_dis = self.n_dis + 1
    #     if self.n_dis > 3:
    #         self.tw_len = {}
    #         self.J = {}
    #         self.n_dis = 0
        
    #     new_window = tk.Toplevel(self.root)
    #     new_window.title('Discriminability convolution')
    #     plot_dis = pt.Plottool(new_window, self.tw_len, self.J_conv)
        
    def discrimin_nd_shift(self):
        timestep = 0.05
        len_t = round(1.2/0.05)
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
                idx0 = 0
                while idx0 + len_t < len(t):
                    E1_w = E1[idx0:idx0 + len_t]
                    E2_w = E2[idx0:idx0 + len_t]
                    # peak1 = max(abs(E1_w))
                    J_temp = 2*simps(abs(E1_w**2-E2_w**2))/simps(abs(E1_w)**2+abs(E2_w)**2)
                    self.tw_len[newkey].append(t[idx0])
                    self.J[newkey].append(J_temp)
                    self.J_base[newkey].append(simps((abs(E1_w**2-E2_w**2)), dx=timestep)/1.2)
                    self.J_log[newkey].append(np.log10(J_temp))
                    idx0 = idx0 + 1
                self.tw_len[newkey] = np.array(self.tw_len[newkey])
                self.J[newkey] = np.array(self.J[newkey])
                self.J_base[newkey] = np.array(self.J_base[newkey])
                self.J_log[newkey] = np.array(self.J_log[newkey])
            # calculate discriminability of superoscillations
            for j in range(0, len(self.idx_so_0)):
                newkey = 'SO_discriminability_' + str(self.n_dis) + '_' +self.bf.find_str_common(self.filenames[self.idx_so_0[j]], self.filenames[self.idx_so_1[j]])
                self.J[newkey] = []
                self.J_base[newkey] = []
                self.J_log[newkey] = []
                self.tw_len[newkey] = []
                t = self.t[self.filenames[self.idx_so_0[j]]]
                E1 = self.x[self.filenames[self.idx_so_0[j]]]
                E2 = self.x[self.filenames[self.idx_so_1[j]]]
                idx0 = 0
                while idx0 + len_t < len(t):
                    E1_w = E1[idx0:idx0 + len_t]
                    E2_w = E2[idx0:idx0 + len_t]
                    # peak1 = max(abs(E1_w))
                    J_temp = 2*simps(abs(E1_w**2-E2_w**2))/simps(abs(E1_w)**2+abs(E2_w)**2)
                    self.tw_len[newkey].append(t[idx0])
                    self.J[newkey].append(J_temp)
                    self.J_base[newkey].append(simps((abs(E1_w**2-E2_w**2)), dx=timestep)/1.2)
                    self.J_log[newkey].append(np.log10(J_temp))
                    idx0 = idx0 + 1
                self.tw_len[newkey] = np.array(self.tw_len[newkey])
                self.J[newkey] = np.array(self.J[newkey])
                self.J_base[newkey] = np.array(self.J_base[newkey])
                self.J_log[newkey] = np.array(self.J_log[newkey])
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
            
        
                
    def Local_f(self):
        newwindow = tk.Toplevel()
        if len(self.idx_so_0) != 0:
            for i in self.idx_so_0:
                key = self.filenames[i]
                newkey = 'local f'+key
                E_cx = hilbert(self.x[key])
                t = self.t[key]
                phase = np.unwrap(np.angle(E_cx))
                self.x[newkey] = self.bf.array_getderivative(t, phase)
                self.t[newkey] = t
        if len(self.idx_so_1) != 0:
            for j in self.idx_so_1:
                key = self.filenames[j]
                E_cx = hilbert(self.x[key])
                t = self.t[key]
                phase = np.unwrap(np.angle(E_cx))
                self.x[newkey] = self.bf.array_getderivative(t, phase)
                self.t[newkey] = t
        plot_localf = pt.Plottool(newwindow, self.t, self.x)
    
    def display_timedelays(self):
        self.td1.set(str(self.so_phase_go.all_time_delays[self.n0.get()]))
        self.td2.set(str(self.so_phase_go.all_time_delays[self.n0.get()+1]))
        self.td3.set(str(self.so_phase_go.all_time_delays[self.n0.get()+2]))
         
    def calculate_phase(self):
        self.ph0 = np.array([self.phase_05.get(),self.phase_06.get(),self.phase_07.get(),self.phase_08.get()])
        self.ph1.set(str(np.array(self.so_phase_go.all_time_delays[self.n0.get()])+self.ph0))
        self.ph2.set(str(np.array(self.so_phase_go.all_time_delays[self.n0.get()+1])+self.ph0))
        self.ph3.set(str(np.array(self.so_phase_go.all_time_delays[self.n0.get()+2])+self.ph0)) 
        
    def plot_all(self):
        newwindow = tk.Toplevel()
        self.so_plotgo = pt.Plottool(newwindow, self.t, self.x)
        
    
    def plot_predictions(self):
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
            plot0.errorbar(
                long_time_window, 
                self.so_phase_go.get_combined_field(time_delays, long_time_window), 
                yerr = self.so_phase_go.get_err_combined_field(time_delays, long_time_window),
                label="theory {}".format(n)
            )

    
        plot0.set_xlabel('time (ps)')
        plot0.set_ylabel('filed (arb. units)')
        plot0.legend()

    
        for n, time_delays in enumerate(all_time_delays[:3]):
            line = plot1.plot(
                self.so_phase_go.time_window, self.so_phase_go.get_combined_field(time_delays, self.so_phase_go.time_window), label="theory {}".format(n))
            line = line[0]
    
            # display raw points
            plot1.errorbar(
                self.so_phase_go.time_window_raw, 
                self.so_phase_go.get_combined_field(time_delays, self.so_phase_go.time_window_raw), 
                yerr = self.so_phase_go.get_err_combined_field(time_delays, self.so_phase_go.time_window_raw),
                color = line.get_color(),
                marker='*',
                linestyle=' ',
            )
    
        plot1.twinx().plot(
            self.so_phase_go.time_window, 
            self.so_phase_go.pulses[self.so_phase_go.largest_freq].interp_field(self.so_phase_go.time_window), 
            'k-',
            label=self.so_phase_go.largest_freq,
        )
       
        plot1.set_xlabel('time (ps)')
        plot1.set_ylabel('filed (arb. units)')
        plot1.legend()
    
        self.canvas0.figure.set_tight_layout(True)
        self.canvas0.draw()
        



        
    
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
    def __init__(self, path, filenames_wf_base):
        self.pulses = {}
        self.max_ampl = []
        # make input variable into dict format
        
        #load data
        for name in filenames_wf_base:
            # load all pulses for the given frequency
            self.filenames = glob.glob(path + name + '_*[0-9].dat')
            self.all_data = [np.loadtxt(_).T[:-1] for _ in self.filenames]
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
        
        #generate initial guesses
        np.random.seed(3112022)

        self.bounds = [(_.peaks_time.min(), _.peaks_time.max()) for _ in self.pulses.values()]
        
        #bounds_mid = [(_.peaks_time.min()+ _.peaks_time.max())/2 for _ in pulses.values()]
        #print(bounds_mid)
        self.rand_initial_guess = np.array([
            np.random.uniform(_[0], _[1], 200 * cpu_count()) for _ in self.bounds
        ])
        self.rand_initial_guess = set(tuple(_) for _ in self.rand_initial_guess.T)
        #print(rand_initial_guess)
        # Pool()
        pool = Pool(processes = int(cpu_count()/2))
        
        self.gradient_descent_results = set(pool.map(self.local_minimization, self.rand_initial_guess))
            
        # self.gradient_descent_results = set(self.local_minimization(_) for _ in (self.rand_initial_guess))
        
       
        self.gradient_descent_results = sorted(self.gradient_descent_results)
        self.intensity_without_ampl_modulation, self.all_time_delays = zip(*self.gradient_descent_results)
    
    def get_combined_field(self,time_delays, time_window):
        return sum(_.interp_field(time_window - delay) for delay, _ in zip(time_delays, self.pulses.values()))
    
    

    def inegral_without_ampl_modulation(self,time_delays):
        return simps(self.get_combined_field(time_delays, self.time_window) ** 2, dx=self.dx)
    
    def integral_inverse(self, time_delays):
        return 1/simps(self.get_combined_field(time_delays, self.time_window) ** 2, dx=self.dx)

    def jac_inegral_without_ampl_modulation(self,time_delays):
        
        E = self.get_combined_field(self,time_delays, self.time_window)
    
        return np.array([     
        -2. * simps(E * _.interp_field.derivative()(self.time_window - delay), dx=self.dx) 
        for delay, _ in zip(time_delays, self.pulses.values())
    ])
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
        return self.inegral_without_ampl_modulation(time_delays), tuple(time_delays)

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
   
    def get_unique_filename(fname):
        return fname.format(datetime.now().strftime("_%m-%d_%H-%M"))    
    
    

        

