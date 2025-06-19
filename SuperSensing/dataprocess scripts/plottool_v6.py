# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 11:03:52 2023

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt
import copy
import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
import tkinter as tk
from tkinter import ttk
import scipy as sp
import pandas as pd
import dataprocessor as dp
from dataprocessor import Basic_functions as bf
from scipy.optimize import curve_fit
import scipy.signal as ss
import curfittool as cft
import configparser as cfp


class Plottool:
    def __init__(self,frame,xvalues,yvalues, nrow = 1, ncol = 1, sort_var = None):
        '''
        

        Parameters
        ----------
        root : tk.Tk()
            DESCRIPTION: tk root to initiate gui base window
        xvalues : dict or list of 1d array or ndarray
            DESCRIPTION: x values to be plotted
        yvalues : dict or list of 1d array or ndarray
            DESCRIPTION: y values to be plotted

        Returns
        -------
        plot of desired axis limit, labels and fonts

        '''
        
        '''global variables'''
        
        self.x = self.formatinput(xvalues)
        self.y = self.formatinput(yvalues)
        self.x_plot = {}
        self.y_plot = {}
        self.frame = frame
        n_color = 0
        for key in list(self.y):
            # self.x_plot[key] = {'data': self.x[key],
            #                     'linestyle': 'solid',
            #                     'linecolor': None,
            #                     'label': key}
            
            self.y_plot[key] = {'linestyle': 'solid',
                                'linecolor': 'C'+str(n_color),
                                'scale': 1,
                                'label': key}
            n_color = n_color+1
        self.action = tk.StringVar(value = 'ready')

                    
        try:
            open('plottool_v5.ini','r')
        except:
            self.autoxlim = tk.IntVar(value = 1)
            self.autoylim = tk.IntVar(value = 1)
            self.fontsize_ax = tk.DoubleVar(value = 20)
            self.fontsize_lg = tk.DoubleVar(value = 20)
            self.labelx = tk.StringVar()
            self.labely = tk.StringVar()
            self.xlim_l = tk.DoubleVar()
            self.xlim_h = tk.DoubleVar()
            self.ylim_l = tk.DoubleVar()
            self.ylim_h = tk.DoubleVar()
            self.legendnames = tk.StringVar()
            self.linestyles = tk.StringVar()
            self.usercolor_varibable = tk.StringVar()
        else:
            pt_config = cfp.ConfigParser()
            pt_config.read('plottool_v5.ini')
            self.autoxlim = tk.IntVar(value = 1)
            self.autoylim = tk.IntVar(value = 1)
            self.fontsize_ax = tk.DoubleVar(value = float(pt_config['init']['fontsize_ax']))
            self.fontsize_lg = tk.DoubleVar(value = float(pt_config['init']['fontsize_lg']))
            self.labelx = tk.StringVar(value = pt_config['init']['label_x'])
            self.labely = tk.StringVar(value = pt_config['init']['label_y'])
            self.xlim_l = tk.DoubleVar(value = float(pt_config['init']['xlim_l']))
            self.xlim_h = tk.DoubleVar(value = float(pt_config['init']['xlim_h']))
            self.ylim_l = tk.DoubleVar(value = float(pt_config['init']['ylim_l']))
            self.ylim_h = tk.DoubleVar(value = float(pt_config['init']['ylim_h']))
            self.legendnames = tk.StringVar(value = pt_config['init']['legendnames'])
            self.linestyles = tk.StringVar(value = pt_config['init']['linestyles'])
            self.usercolor_varibable = tk.StringVar(value = pt_config['init']['colors'])
            
        self.nrow = nrow
        self.ncol = ncol
        self.nrow_new = tk.IntVar(value = 1)
        self.ncol_new = tk.IntVar(value = 1)
        self.plot_num = tk.IntVar()
        self.figsize = tk.StringVar(value = 'medium')
        self.all_plot_nums = np.arange(0, nrow*ncol, step = 1).tolist()
        self.idx_selectx = tuple([])
        self.idx_selecty = tuple([])
        self.auto_selecty = tk.IntVar(value = 1)
        self.line_weight = tk.IntVar(value = 1)
        self.config_largeplot = tk.StringVar()
        self.fig_name = tk.StringVar(value = 'Fig.1')
       
        self.user_linestyle_switch = tk.IntVar(value = 0)
        self.legendnames_display = tk.StringVar(value = 'None')
        self.lg_pos = tk.StringVar(value = 'best')
        self.userlgname = tk.IntVar(value=0)
        self.plot_scale = tk.StringVar()
        self.plot_scale_switch = tk.IntVar(value = 0)
        self.sort_keyword = tk.StringVar(value = 'K')
        
        self.names_y = sorted(list(self.y), reverse = True)
        self.names_x = sorted(list(self.x), reverse = True)   
            
        self.listvar_x = tk.StringVar(value=self.names_x)
        self.listvar_y = tk.StringVar(value=self.names_y)
        self.text_weight = tk.StringVar(value = 'normal')
        
        self.legendon = tk.IntVar(value=1)
        self.usercolor_switch = tk.IntVar(value=0)
        self.usercolor_value = []
        self.gridswitch = tk.IntVar(value = 1)
        self.color_type = tk.StringVar(value = 'name')
        
        self.selection_x = tk.StringVar(value='Selected x spectrum: \n None')
        self.selection_y = tk.StringVar(value='Selected y spectrum: \n None')
        
        
        '''first window for showing figure'''
        
        self.window1 = ttk.Labelframe(self.frame,text='plot')
        self.window1.grid(column=0,row=0, columnspan = 3, sticky='w')
        self.figure, self.ax = plt.subplots(self.nrow, self.ncol, figsize = (10, 5), dpi = 100)
        self.figure.set_tight_layout(True)
        self.canvas0 = FigureCanvasTkAgg(figure = self.figure, master=self.window1)
        # self.canvas0.figure, self.ax = plt.subplots(self.nrow, self.ncol, figsize = (15, 8))
        self.canvas0.get_tk_widget().grid(column = 0,row = 0,columnspan = 16,rowspan = 10, sticky='w')
        self.canvas0.draw()
        button_operate_plot = ttk.Button(self.window1,text='examine',command=self.exam_plot)
        button_curvefit = ttk.Button(self.window1, text = 'curve fit', command = self.start_curvefit)
        self.button_plot = ttk.Button(self.window1,text='update plot',command=self.plotfig)
        
        label_fig_name = ttk.Label(self.window1, text = 'Enter Figure name here: ')
        entry_fig_name = ttk.Entry(self.window1, textvariable = self.fig_name)
        button_save_figinfo= ttk.Button(self.window1, text = 'Save fig info', command = self.save_fig_info)
      
        
        button_exportdata = ttk.Button(self.window1,text='Export Data as shown in plot',command=self.export_data)
        entry_config_largeplot = ttk.Entry(self.window1, textvariable = self.config_largeplot, width = 30)
        button_large_plot = ttk.Button(self.window1, text = 'make large plot', command = self.make_large_plot)
        
        
        
        label_combobox_nrow = ttk.Label(self.window1, text = 'rows of subplots: ')
        label_combobox_ncol = ttk.Label(self.window1, text = 'columns of subplots: ')
        self.combobox_nrow = ttk.Combobox(self.window1, textvariable = self.nrow_new, width = 5)
        self.combobox_nrow['values'] = [1, 2, 3]
        self.combobox_nrow.bind('<<ComboboxSelected>>', self.create_subplots)
        self.combobox_ncol = ttk.Combobox(self.window1, textvariable = self.ncol_new, width = 5)
        self.combobox_ncol['values'] = [1, 2, 3]
        self.combobox_ncol.bind('<<ComboboxSelected>>', self.create_subplots)
        label_combobox_plotnum = ttk.Label(self.window1, text = 'Select subplot number to plot: ', width = 40)
        self.combobox_plot_num = ttk.Combobox(self.window1, textvariable = self.plot_num, width = 5)
        self.combobox_plot_num['values'] = self.all_plot_nums
        
        
        button_operate_plot.grid(column = 0,row = 12, sticky = 'w')
        button_curvefit.grid(column = 1, row = 12, sticky = 'w')
        button_exportdata.grid(column = 2,row = 12, columnspan = 2, sticky = 'w')
        entry_config_largeplot.grid(column = 4, row = 12, columnspan = 2, sticky = 'w')
        button_large_plot.grid(column = 6, row = 12, sticky = 'w')
        label_combobox_nrow.grid(column = 0, row = 13, sticky = 'w')
        self.combobox_nrow.grid(column = 1, row = 13, sticky = 'w')
        label_combobox_ncol.grid(column = 2, row = 13, sticky = 'w')
        self.combobox_ncol.grid(column = 3, row = 13, sticky = 'w')
        label_combobox_plotnum.grid(column = 4, row = 13, columnspan = 4, sticky = 'w')
        self.combobox_plot_num.grid(column = 8, row = 13, sticky = 'w')
        
        self.button_plot.grid(column = 9, row = 13, sticky = 'w')
        
        label_fig_name.grid(column = 0, row = 14, sticky = 'w')
        entry_fig_name.grid(column = 1, row = 14, sticky = 'w')
        button_save_figinfo.grid(column = 2, row = 14, sticky = 'w')
        
        '''second window for adjusting figure'''
    
        self.window = ttk.Labelframe(frame,text='plot editor')
        self.window.grid(column=0,row=1, columnspan = 2, sticky = 'w')
        
        
        self.entry_xlabel = ttk.Entry(self.window,textvariable=self.labelx,width=30)
        self.entry_ylabel = ttk.Entry(self.window,textvariable=self.labely,width=30)
        self.entry_fontsize_lg = ttk.Entry(self.window,textvariable=self.fontsize_lg,width=10)
        self.entry_fontsize_ax = ttk.Entry(self.window,textvariable=self.fontsize_ax,width=10)
        label_line_weight = ttk.Label(self.window, text = '  Line width: ')
        entry_line_weight = ttk.Entry(self.window, textvariable = self.line_weight, width = 10)
        
        
        label_xlim1 = ttk.Label(self.window,text='x from: ')
        label_xlim2 = ttk.Label(self.window,text='to: ')
        label_ylim1 = ttk.Label(self.window,text='y from: ')
        label_ylim2 = ttk.Label(self.window,text='to: ')
        entry_xl = ttk.Entry(self.window,textvariable=self.xlim_l,width=10)
        entry_xh = ttk.Entry(self.window,textvariable=self.xlim_h,width=10)
        entry_yl = ttk.Entry(self.window,textvariable=self.ylim_l,width=10)
        entry_yh = ttk.Entry(self.window,textvariable=self.ylim_h,width=10)
        
        checkbutton_autoxlim = ttk.Checkbutton(self.window,variable=self.autoxlim, text = 'auto x limits')
        checkbutton_autoylim = ttk.Checkbutton(self.window,variable=self.autoylim, text = 'auto y limits')
        checkbutton_legendon = ttk.Checkbutton(self.window,variable=self.legendon, text = 'legend on  ')
        
        label_combobox_textweight = ttk.Label(self.window, text = '  text weight:')
        self.combobox_textweight = ttk.Combobox(self.window, textvariable = self.text_weight, width = 10)
        self.combobox_textweight['value'] = ('bold', 'light', 'heavy', 'normal')
        self.combobox_textweight.bind('<<ComboboxSelected>>', self.new_parem_plot)
        
        label_xlabel = ttk.Label(self.window,text='x-axis label: ')
        label_ylabel = ttk.Label(self.window,text='y-axis label: ')
        self.label_font_lg = ttk.Label(self.window,text='  legend font size: ')
        self.label_font_ax = ttk.Label(self.window,text='  axis font size: ')
        # self.label_grid = ttk.Label(self.window,text='grid on: ')
        self.checkbutton_grid = ttk.Checkbutton(self.window,variable=self.gridswitch, text = 'grid on  ')
        
        # widges for legends
        label_entry_lgnames = tk.Label(self.window, text = 'Enter legends here, use + as seperator:')
        self.entry_legendnames = ttk.Entry(self.window,textvariable=self.legendnames,width = 70)
        # label_lg_0 = ttk.Label(self.window4, text = 'Following legends will be used: ')
        # label_show_loaded_lgs = ttk.Label(self.window4, textvariable = self.legendnames_display)
        # label_lg_pos = ttk.Label(self.window4, text = 'Legend position: ')
        combobox_lg_pos  = ttk.Combobox(self.window, textvariable = self.lg_pos, width = 10)
        combobox_lg_pos['value'] = ['best', 'upper right', 'lower right', 'upper left', 'lower left', 'center left', 'center right',
                                    'upper center', 'lower center']
        button_loadlegendnames = ttk.Button(self.window,text='Load',command=self.load_lg_names)
        
        
        
        
        label_xlabel.grid(column = 0, row = 0, columnspan = 1, sticky = 'w')
        self.entry_xlabel.grid(column = 1, row = 0, columnspan = 3, sticky = 'w')
        label_xlim1.grid(column = 4,row = 0, sticky = 'w')
        entry_xl.grid(column = 5,row = 0, sticky = 'w')  
        label_xlim2.grid(column = 6,row = 0, sticky = 'w')
        entry_xh.grid(column = 7,row = 0, sticky = 'w') 
        checkbutton_autoxlim.grid(column = 8,row = 0, columnspan = 2, sticky = 'w')
        label_ylabel.grid(column = 0, row = 1, sticky = 'w')
        self.entry_ylabel.grid(column = 1, row = 1, columnspan = 3, sticky = 'w')
        label_ylim1.grid(column = 4,row = 1, sticky = 'w')
        entry_yl.grid(column = 5,row = 1, sticky = 'w')
        label_ylim2.grid(column = 6,row = 1, sticky = 'w')
        entry_yh.grid(column = 7,row = 1, sticky = 'w')
        checkbutton_autoylim.grid(column = 8,row = 1, columnspan = 2, sticky = 'w')
        
        
        self.label_font_lg.grid(column = 0, row = 2, columnspan = 1, sticky = 'w')
        self.entry_fontsize_lg.grid(column = 1,row = 2, sticky = 'w')
        self.label_font_ax.grid(column = 2, row = 2, columnspan = 1, sticky = 'w')
        self.entry_fontsize_ax.grid(column = 3,row = 2, sticky = 'w')
        label_line_weight.grid(column = 4, row = 2, sticky = 'w')
        entry_line_weight.grid(column = 5, row = 2, sticky = 'w')
        self.checkbutton_grid.grid(column = 6,row = 2, sticky = 'w')
        checkbutton_legendon.grid(column = 7,row = 2, sticky = 'w')
        label_combobox_textweight.grid(column  = 8, row = 2, sticky = 'w')
        self.combobox_textweight.grid(column = 9, row = 2, sticky = 'w')
        
        label_entry_lgnames.grid(column = 0, row = 3, columnspan = 2, sticky = 'w')
        self.entry_legendnames.grid(column = 0, row = 4, columnspan = 5, sticky = 'w')
        # label_lg_pos.grid(column = 4, row = 1, columnspan = 1, sticky = 'w')
        combobox_lg_pos.grid(column = 6, row = 4, columnspan = 1, sticky = 'w')
        button_loadlegendnames.grid(column = 7, row = 4, columnspan = 1, sticky = 'w')
        
        
        '''third window for choose x and y to plot'''
        
        
        # main window
        self.window3 = ttk.Labelframe(self.frame, text = 'select x and y')
        self.window3.grid(column = 3, row = 0, rowspan = 2, sticky = 'nw')
        
        label_show_select_x = ttk.Label(self.window3,textvariable=self.selection_x)
        label_show_select_y = ttk.Label(self.window3,textvariable=self.selection_y)
        self.ydata_select = tk.Listbox(self.window3,listvariable=self.listvar_y,selectmode='extended',width=60, state = 'disabled')
        self.xdata_select = tk.Listbox(self.window3,listvariable=self.listvar_x,selectmode='extended',width=60)
        button_loadselection_x = ttk.Button(self.window3,text='load x selection',command=self.loadselection_x)
        button_loadselection_y = ttk.Button(self.window3,text='load y selection',command=self.loadselection_y)
        checkbutton_auto_selecty = ttk.Checkbutton(self.window3, text = 'auto select y on: ', 
                                                   variable = self.auto_selecty, command = self.config_listbox_ydata)
        
        button_sortlist = ttk.Button(self.window3,text='sort x and y lists', command = self.sort_xy_list)
        cbox_sort_keyword = ttk.Combobox(self.window3, textvariable = self.sort_keyword, width = 10)
        cbox_sort_keyword['value'] = ['K', 'T', 'D']
        
        cbox_sort_keyword.grid(column = 0, row = 0, sticky = 'w')
        button_sortlist.grid(column = 1, row = 0, sticky = 'w')
        self.xdata_select.grid(column = 0, row = 1, columnspan = 2, sticky = 'w')
        button_loadselection_x.grid(column = 0, row = 2, sticky = 'w')
        self.ydata_select.grid(column = 0, row = 3, columnspan = 2, sticky = 'w')
        checkbutton_auto_selecty.grid(column = 0, row = 4, sticky = 'w')
        button_loadselection_y.grid(column = 0, row = 5, sticky = 'w')
        label_show_select_x.grid(column  = 0, row = 6, columnspan = 2, sticky = 'w')
        label_show_select_y.grid(column  = 1, row = 6, columnspan = 2, sticky = 'w')
        
        
        
        '''forth window to change detail of the figure'''
        
        
        self.window4 = ttk.Labelframe(self.frame, text = 'Edit Curves')
        self.window4.grid(column = 2, row = 1, columnspan = 2, sticky = 'nw')
        
        
        
    
        # widges for colors
        label_entry_colors = tk.Label(self.window4, text = 'Enter colors here, use + as seperator:')
        button_loadcolors = ttk.Button(self.window4,text='load colors',command = self.load_colors)
        self.usercolor_entry = ttk.Entry(self.window4, textvariable = self.usercolor_varibable, width = 70)
        cbox_color_type = ttk.Combobox(self.window4, textvariable = self.color_type, width = 10)
        cbox_color_type['value'] = ['name', 'RGB']
        
        # widges for curve scales
        label_entry_scales = tk.Label(self.window4, text = 'Enter curve scales here, use + as seperator:')
        entry_plot_scales = ttk.Entry(self.window4, textvariable = self.plot_scale, width = 70)
        button_plotscales_load = ttk.Button(self.window4, text = 'load', command = self.load_plot_scale, width = 10)
        
        # widges for linestyle
        label_entry_linestyle = tk.Label(self.window4, text = 'Enter linestyles here, use + as seperator:')
        entry_linestyles = ttk.Entry(self.window4, textvariable = self.linestyles, width = 70)
        button_user_linestyle = ttk.Button(self.window4, text = 'load', command = self.load_linestyles, width = 10)
        
        
        button_config = ttk.Button(self.window4, text = 'save as default', command = self.build_config)
        
        
        # grid widges
        
        
        
        label_entry_colors.grid(column = 0, row = 2, columnspan = 2, sticky = 'w')
        self.usercolor_entry.grid(column = 0, row = 3, columnspan = 4, sticky = 'w')
        cbox_color_type.grid(column = 4, row = 3, sticky = 'w')
        button_loadcolors.grid(column = 5, row = 3, sticky = 'w')
        
        label_entry_scales.grid(column = 0, row = 4, columnspan = 2, sticky = 'w')
        entry_plot_scales.grid(column = 0, row = 5, columnspan = 4, sticky = 'w')
        button_plotscales_load.grid(column = 4, row = 5, columnspan = 1, sticky = 'w')
        
        label_entry_linestyle.grid(column = 0, row = 6, columnspan = 2, sticky = 'w')
        entry_linestyles.grid(column = 0, row = 7, columnspan = 4, sticky = 'w')
        button_user_linestyle.grid(column = 4, row = 7, columnspan = 1, sticky = 'w')
        
        button_config.grid(column = 0, row = 8, columnspan = 1, sticky = 'w')
        
        
        '''
        bottom window for status bar
        '''
        self.window6 = tk.Frame(frame)
        self.window6.grid(column=0,row = 2, sticky = 'w')

        label_status = ttk.Label(self.window6,text='status: ', font = ('Times', 15))
        label_showaction = ttk.Label(self.window6,textvariable=self.action, font=('Times', 15), width = 30)
        
        label_status.grid(column=0,row=0, sticky = 'w')
        label_showaction.grid(column=1,row=0, sticky = 'w')
  
        
    def build_config(self):
        pt_config = cfp.ConfigParser()
        pt_config['init'] = {'label_x': self.labelx.get(),
                             'label_y': self.labely.get(),
                             'fontsize_ax': str(self.fontsize_ax.get()),
                             'fontsize_lg': str(self.fontsize_lg.get()),
                             'xlim_l': str(self.xlim_l.get()),
                             'xlim_h': str(self.xlim_h.get()),
                             'ylim_l': str(self.ylim_l.get()),
                             'ylim_h': str(self.ylim_h.get()),
                             'legendnames': self.legendnames.get(),
                             'linestyles': self.linestyles.get(),
                             'colors': self.usercolor_varibable.get()}
        with open('plottool_v5.ini','w') as plottool_v5_init:
            pt_config.write(plottool_v5_init)
        self.action.set('current entry values are saved as default values')
    

    
    def loadselection_x(self):
        self.selection_x.set('Selected x spectrums: \n')
        self.idx_selectx = self.xdata_select.curselection()
        for i in self.idx_selectx:
            self.selection_x.set(self.selection_x.get()+self.names_x[i]+'\n')
        if self.auto_selecty.get() == 1:
            self.idx_selecty = self.idx_selectx
            self.selection_y.set('Selected y spectrums: \n')
            for j in self.idx_selecty:
                self.selection_y.set(self.selection_y.get()+self.names_y[j]+'\n')
            self.action.set('x and y data selected') 
    def loadselection_y(self):
        self.idx_selecty = self.ydata_select.curselection()
        self.selection_y.set('Selected y spectrums: \n')
        for i in self.idx_selecty:
            self.selection_y.set(self.selection_y.get()+self.names_y[i]+'\n')
            
        self.action.set('y data selected')
        
    def sort_xy_list(self):
        if self.sort_keyword.get() == 'K':
            self.names_x.sort(key = bf.dict_key_get_T)
            self.names_y.sort(key = bf.dict_key_get_T)
            self.listvar_x.set(self.names_x)
            self.listvar_y.set(self.names_y)
        if self.sort_keyword.get() == 'T':
            self.names_x.sort(key = bf.dict_key_get_B)
            self.names_y.sort(key = bf.dict_key_get_B)
            self.listvar_x.set(self.names_x)
            self.listvar_y.set(self.names_y)
        if self.sort_keyword.get() == 'D':
            self.names_x.sort(key = bf.dict_key_get_D)
            self.names_y.sort(key = bf.dict_key_get_D)
            self.listvar_x.set(self.names_x)
            self.listvar_y.set(self.names_y)
            
            
    def exam_plot(self):
        root = tk.Toplevel()
        frame = ttk.Frame(root)
        frame.pack(fill=tk.BOTH,expand=True)
        self.canvas1 = FigureCanvasTkAgg(self.canvas0.figure,master=frame)
        self.canvas1.figure = self.canvas0.figure
        self.canvas1.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.canvas1.draw()
        self.toolbar = NavigationToolbar2Tk(self.canvas1,frame)
        
        
    def start_curvefit(self):
        new_root = tk.Toplevel(self.frame)
        len_data = len(self.idx_selectx)
        if len_data == 1:
            x = self.x[self.names_x[self.idx_selectx[0]]]
            y = self.y[self.names_y[self.idx_selecty[0]]]
            idx0 = np.where(x>self.xlim_l.get())[0][0]
            idx1 = np.where(x<self.xlim_h.get())[0][-1]
            x = x[idx0:idx1]
            y = y[idx0:idx1]
            self.curfitgo = cft.Curvefit(new_root, x, y)
        else:
            self.action.set('Please only select 1 pair of data to do curve fitting')
            
            
    def config_listbox_ydata(self):
        if self.auto_selecty.get() == 1:
            self.ydata_select.config(state = 'disabled')
        else:
            self.ydata_select.config(state = 'normal')
             
             
    def formatinput(self,x):
        '''
        

        Parameters
        ----------
        x : list or dict or ndarray
            input data to be format into a dict type

        Returns
        -------
        x_new : dict
            output dict type data

        '''
        if type(x) is list:
            x_new = dict()
            for i in range(len(x)):
                name = 'x'+str(i)
                x_new[name] = x[i]
                
        if type(x) is np.ndarray:
            x_new = dict()
            if np.size(x) == len(x):
                    name = 'x0'
                    x_new[name] = x
            else:
                for i in range(np.size(x,1)):
                    name = 'x'+str(i)
                    x_new[name] = x[:,i]
                    
                    
        if type(x) is dict:
            x_new = copy.deepcopy(x)
        
        if type(x) is pd.core.frame.DataFrame:
            x_new = x
        
        return x_new
    
    
    def new_parem_plot(self, comboboxx_event):
        mpl.rc('font', weight = self.text_weight.get())
        
    
    def create_subplots(self, combobox_event):
        self.ncol = self.ncol_new.get()
        self.nrow = self.nrow_new.get()
        self.all_plot_nums = (np.arange(0, self.nrow*self.ncol, step = 1).tolist())
        self.combobox_plot_num['values'] = self.all_plot_nums
        plt.close()
        self.figure, self.ax = plt.subplots(self.nrow, self.ncol, figsize = (18, 10))
        self.canvas0.figure = self.figure
        self.canvas0.draw()
        
    def load_lg_names(self):
        self.userlegendnames = self.legendnames.get().split('+')
        self.legendnames_display.set(self.legendnames.get().replace('+', '\n'))
        if len(self.userlegendnames) == len(self.idx_selecty):
            for i, j in enumerate(self.idx_selecty):
                key = self.names_y[j]
                self.y_plot[key]['label'] = self.userlegendnames[i]
            self.action.set('user defined legend names loaded')
        else:
            self.action.set('Error! number of legend names entered does not match number of selected x contents')
            
    def load_colors(self):
        if self.color_type.get() == 'name':
            self.loadcolors_name()
        if self.color_type.get() == 'RGB':
            self.loadcolors_tuple()
    
    
    def loadcolors_name(self):
        self.usercolor_value = self.usercolor_varibable.get().split('+')
        if len(self.usercolor_value) == len(self.idx_selecty):
            for i,j in enumerate(self.idx_selecty):
                key = self.names_y[j]
                self.y_plot[key]['linecolor'] = self.usercolor_value[i]
            self.action.set('colors loaded')
        else:
            self.action.set('Error! number of colors entered does not match number of selected x contents')
        
    def loadcolors_tuple(self):
        self.usercolor_value = self.usercolor_varibable.get().split('+')
        for i in range(0,len(self.usercolor_value)):
            self.usercolor_value[i] = eval(self.usercolor_value[i])
        if len(self.usercolor_value) == len(self.idx_selecty):
            for i in range(len(self.idx_selecty)):
                key = self.names_y[i]
                self.y_plot[key]['linecolor'] = self.usercolor_value[i]
            self.action.set('colors loaded')
        else:
            self.action.set('Error! number of colors entered does not match number of selected y contents')
        self.action.set('colors loaded')
        
    def load_plot_scale(self):
        self.plot_scale_values = self.plot_scale.get().split('+')
        if len(self.plot_scale_values) == len(self.idx_selecty):
            for i,j in enumerate(self.idx_selecty):
                key = self.names_y[j]
                self.y_plot[key]['scale'] = float(self.plot_scale_values[i])
            self.action.set('plot scales loaded')
        else:
            self.action.set('Error! number of plot scale values entered does not match number of selected y contents')
            
            
    def load_linestyles(self):
        self.linestyles_load = self.linestyles.get().split('+')
        if len(self.linestyles_load) == len(self.idx_selecty):
            for i,j in enumerate(self.idx_selecty):
                key = self.names_y[j]
                self.y_plot[key]['linestyle'] = self.linestyles_load[i]
            self.action.set('linestyles loaded')
        else:
            self.action.set('Error! number of linestyles entered does not match number of selected y contents') 
            
   
            
    def plotfig(self):
        # make sure x data and y data are selected
        if len(self.idx_selectx) == 0 or len(self.idx_selecty) == 0:
            self.action.set('Error! please select x and y to be plotted')
            # self.button_plot.state(['disabled'])
            return
        
        if type(self.ax) is np.ndarray:
            if len(np.shape(self.ax)) == 1: 
                plot0 = self.ax[self.plot_num.get()]
            else:
                nrow = int(np.floor(self.plot_num.get()/self.ncol))
                ncol = int(self.plot_num.get()%self.ncol)
                plot0 = self.ax[nrow, ncol]
        else:
            plot0 = self.ax
       
        plot0.cla()
        # if self.usercolor_switch.get() == 1:
        if len(self.idx_selectx) == 1:
            key_x = self.names_x[self.idx_selectx[0]]
            for i in self.idx_selecty:
                key_y = self.names_y[i]
                scale = self.y_plot[key_y]['scale']
                if len(self.x[key_x]) >= len(self.y[key_y]):
                    len_plot = len(self.y[key_y])
                else:
                    len_plot = len(self.x[key_x])
                    
                plot0.plot(self.x[key_x][0:len_plot], self.y[key_y][0:len_plot]*scale,
                           label = self.y_plot[key_y]['label'],
                           color = self.y_plot[key_y]['linecolor'], 
                           linestyle = self.y_plot[key_y]['linestyle'],
                           linewidth = self.line_weight.get())
                # plt.tight_layout()
        else:
            if len(self.idx_selectx) == len(self.idx_selecty):
                for i in range(0, len(self.idx_selectx)):
                    key_x = self.names_x[self.idx_selectx[i]]
                    key_y = self.names_y[self.idx_selecty[i]]
                    scale = self.y_plot[key_y]['scale']
                    plot0.plot(self.x[key_x], self.y[key_y]*scale,
                               label = self.y_plot[key_y]['label'], 
                               color = self.y_plot[key_y]['linecolor'],
                               linestyle = self.y_plot[key_y]['linestyle'],
                               linewidth =self.line_weight.get())
                    # plt.tight_layout()
            else:
                self.action.set('Error! please select same number of x and y')
           
       
        if self.autoxlim.get() == 0:
            plot0.set_xlim(left=float(self.xlim_l.get()),right=float(self.xlim_h.get()))
        if self.autoylim.get() == 0:
            plot0.set_ylim(bottom=float(self.ylim_l.get()),top=float(self.ylim_h.get()))
        plot0.set_xlabel(self.labelx.get(),fontsize = self.fontsize_ax.get())
        plot0.set_ylabel(self.labely.get(),fontsize = self.fontsize_ax.get())
        plot0.tick_params(axis='both', labelsize = self.fontsize_ax.get())
        if self.legendon.get() == 1:
            plot0.legend(loc = self.lg_pos.get(), fontsize = self.fontsize_lg.get())
        plot0.grid(self.gridswitch.get())
        # self.canvas0.figure.savefig(path+'\\'+'Waveforms_fd.pdf')
        self.canvas0.draw()
        
    def save_fig_info(self):
        path = tk.filedialog.askdirectory()
        self.fig_info = {}
        if len(self.idx_selectx) == 1:
            key_x = self.names_x[self.idx_selectx[0]]
            for i in self.idx_selecty:
                key_y = self.names_y[i]
                scale = self.y_plot[key_y]['scale']
                # if len(self.x[key_x]) >= len(self.y[key_y]):
                #     len_plot = len(self.y[key_y])
                # else:
                #     len_plot = len(self.x[key_x])
                curve_label = self.y_plot[key_y]['label']
                self.fig_info[curve_label] = np.stack((self.x[key_x], self.y[key_y]), axis = 1)
            # np.savez(path+'/'+self.fig_name.get()+curve_label, x = self.x[key_x], y = self.y[key_y]*scale)
            np.savez(path+'/'+self.fig_name.get(), **self.fig_info)
        else:
            if len(self.idx_selectx) == len(self.idx_selecty):
                for i in range(0, len(self.idx_selectx)):
                    key_x = self.names_x[self.idx_selectx[i]]
                    key_y = self.names_y[self.idx_selecty[i]]
                    scale = self.y_plot[key_y]['scale']
                    curve_label = self.y_plot[key_y]['label']
                    # self.fig_info[curve_label] = {self.labelx.get(): self.x[key_x],
                    #                     self.labely.get(): self.y[key_y]*scale}
                    self.fig_info[curve_label] = np.stack((self.x[key_x], self.y[key_y]), axis = 1)
                # np.savez(path+'/'+self.fig_name.get()+curve_label, x = self.x[key_x], y = self.y[key_y]*scale)
                np.savez(path+'/'+self.fig_name.get(), **self.fig_info)
            else:
                self.action.set('Error! please select same number of x and y')

    def make_large_plot(self):
        window = tk.Toplevel()
        frame = ttk.Frame(window)
        frame.pack(fill=tk.BOTH,expand=True)
        # button_operate_plot = ttk.Button(window,text='examine',command=self.exam_plot)
        # button_operate_plot.grid(column = 0, row = 0)
        if self.config_largeplot.get() == '':
            figure, ax = plt.subplots(1, 1, figsize = (4, 3), dpi = 300)
        else:
            size_col = float(self.config_largeplot.get().split(',')[0])
            size_row = float(self.config_largeplot.get().split(',')[1])
            size_dpi = float(self.config_largeplot.get().split(',')[2])
            figure, ax = plt.subplots(1, 1, figsize = (size_col, size_row), dpi = size_dpi)
        figure.set_tight_layout(True)
        canvas0 = FigureCanvasTkAgg(figure = figure, master=frame)
        
        canvas0.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        plot0 = ax
        if len(self.idx_selectx) == 1:
            key_x = self.names_x[self.idx_selectx[0]]
            for i in self.idx_selecty:
                key_y = self.names_y[i]
                scale = self.y_plot[key_y]['scale']
                if len(self.x[key_x]) >= len(self.y[key_y]):
                    len_plot = len(self.y[key_y])
                else:
                    len_plot = len(self.x[key_x])
                    
                plot0.plot(self.x[key_x][0:len_plot], self.y[key_y][0:len_plot]*scale,
                           label = self.y_plot[key_y]['label'],
                           color = self.y_plot[key_y]['linecolor'], 
                           linestyle = self.y_plot[key_y]['linestyle'],
                           linewidth = self.line_weight.get())
                # plt.tight_layout()
        else:
            if len(self.idx_selectx) == len(self.idx_selecty):
                for i in range(0, len(self.idx_selectx)):
                    key_x = self.names_x[self.idx_selectx[i]]
                    key_y = self.names_y[self.idx_selecty[i]]
                    scale = self.y_plot[key_y]['scale']
                    plot0.plot(self.x[key_x], self.y[key_y]*scale,
                               label = self.y_plot[key_y]['label'], 
                               color = self.y_plot[key_y]['linecolor'],
                               linestyle = self.y_plot[key_y]['linestyle'],
                               linewidth =self.line_weight.get())
        if self.autoxlim.get() == 0:
            plot0.set_xlim(left=float(self.xlim_l.get()),right=float(self.xlim_h.get()))
        if self.autoylim.get() == 0:
            plot0.set_ylim(bottom=float(self.ylim_l.get()),top=float(self.ylim_h.get()))
        plot0.set_xlabel(self.labelx.get(),fontsize = self.fontsize_ax.get())
        plot0.set_ylabel(self.labely.get(),fontsize = self.fontsize_ax.get())
        plot0.tick_params(axis='both', labelsize = self.fontsize_ax.get())
        if self.legendon.get() == 1:
            plot0.legend(loc = self.lg_pos.get(), fontsize = self.fontsize_lg.get(), framealpha = 0.5)
        plot0.grid(self.gridswitch.get())
        canvas0.draw()
        toolbar = NavigationToolbar2Tk(canvas0,frame)
        
    
    def export_data(self):
        '''
        export processed data
        
        '''
        path = tk.filedialog.askdirectory()
        if len(self.idx_selectx) != 1:
            if len(self.idx_selectx) == len(self.idx_selecty):
                for i in range(0, len(self.idx_selectx)):
                    key_x = self.names_x[self.idx_selectx[i]]
                    key_y = self.names_y[self.idx_selecty[i]]
                    newname_x = 'export_' + key_x
                    newname_y = 'export_' + key_y
                    scale = self.y_plot[key_y]['scale']
                    if self.autoxlim.get() == 0:
                        self.x[newname_x], self.y[newname_y] = dp.Basic_functions().array_chop(self.x[key_x],self.y[key_y], self.xlim_l.get(), self.xlim_h.get())
                    else:
                        self.x[newname_x] = self.x[key_x]
                        self.y[newname_y] = self.y[key_y]
                    column_fill = np.zeros_like(self.y[newname_y])
                    data_output = np.stack((self.x[newname_x], scale*self.y[newname_y], column_fill),axis=1)
                    # datasave = open(path+ '/' + newname_x + '.dat', 'w')
                    np.savetxt(path + '/' + newname_x + '_1.dat', data_output)
                    # datasave.close()
                    del newname_x,newname_y,data_output
            else:
                self.action.set('Error! selection same number of x and y')
                        
        else:
            key_x = self.names_x[self.idx_selectx[0]]
            newname_x = 'export_' + key_x
            for i in range(0, len(self.idx_selecty)):
                key_y = self.names_y[self.idx_selecty[i]]
                newname_y = 'export_' + key_y
                scale = self.y_plot[key_y]['scale']
                if self.autoxlim.get() == 0:
                    self.x[newname_x],self.y[newname_y] = dp.Basic_functions().array_chop(self.x[key_x],self.y[key_y], self.xlim_l.get(), self.xlim_h.get())
                else:
                    self.x[newname_x] = self.x[key_x]
                    self.y[newname_y] = self.y[key_y]
                column_fill = np.zeros_like(self.y[newname_y])
                data_output = np.stack((self.x[newname_x], scale*self.y[newname_y], column_fill),axis=1)
                # datasave = open(path + '/'+ newname_y + '.dat', 'w')
                np.savetxt(path + '/' + newname_y + '_1.dat', data_output)
                # datasave.close()
                del newname_x,newname_y, data_output
        self.action.set('data exported')
        

if __name__ == '__main__':
    x = np.linspace(1, 10,10)
    y = x**2
    root=tk.Tk()
    mainframe = ttk.Frame(root)
    mainframe.grid(column = 0, row = 0)
    pt=Plottool(mainframe, x, y)
    root.mainloop()