# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 16:02:24 2023

@author: Admin

"""

import tkinter as tk
from tkinter import ttk
import copy

class Wave_simulator:
    def __init__(self, frame, tavg, xavg):
        self.window = ttk.Labelframe(frame, text = 'Operation Panel')
        self.window.grid(column = 0, row = 0, sticky = 'w')
        self.t = copy.deepcopy(tavg)
        self.x = copy.deepcopy(xavg)
        self.t_backup = copy.deepcopy(tavg)
        self.x_backup = copy.deepcopy(xavg)
        
        self.input_mode = tk.StringVar(value = 'experimental')
        self.filenames = tk.StringVar(value = list(self.x))
        
        label_input_mode = ttk.Label(self.window, text = 'Select input mode: ')
        cbox_input_mode = ttk.Combobox(self.window, textvariable = self.input_mode, width = 10)
        cbox_input_mode['value'] = ['experimental', 'arbitrary']
        
        
        label_input_mode.grid(column = 0, row = 0, sticky = 'w')
        cbox_input_mode.grid(column = 1, row = 0, sticky = 'w')
        
        self.window_exp = ttk.Labelframe(self.window, text = 'Simulate with experimental data')
        self.window_arb = ttk.Labelframe(self.window, text = 'Simulate with arbitrary function')
        
        '''window content of self.window_exp'''
        lbox_filenames = tk.Listbox(self.window_exp, listvariable = self.filenames, selectmode = 'extended')
        
        
        
if __name__ == '__main__':
    root = tk.Tk()
    simugo = Wave_simulator(root)
    root.mainloop()