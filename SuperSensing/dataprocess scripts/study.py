import tkinter as tk
from tkinter import ttk
import dataprocessor as dp
import numpy as np
import plottool_v6 as pt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
import copy
import matplotlib.pyplot as plt
import configparser as cfp
from numpy import pi
from dataprocessor import Basic_functions as bf

class study_transmission:
    def __init__(self,frame,freq,sxsam,sxref):
        self.frame_tr = frame
        self.freq = copy.deepcopy(freq)
        self.sxsam = copy.deepcopy(sxsam)
        self.sxref = copy.deepcopy(sxref)
        self.sam_selected = tk.StringVar(value = 'None \n')
        self.ref_selected = tk.StringVar(value = 'None \n')
        self.window1 = ttk.Labelframe(self.frame_tr,text='main')
        self.window1.grid(column=0,row=0)
        self.fnames_sam = list(sxsam)
        self.fnames_ref = list(sxref)
        self.list_sam = tk.StringVar(value=self.fnames_sam)
        self.list_ref = tk.StringVar(value=self.fnames_ref)
        self.action = tk.StringVar(value='ready')
        self.listbox_sam= tk.Listbox(self.window1,listvariable=self.list_sam,selectmode='extended', width = 30)
        self.listbox_ref= tk.Listbox(self.window1,listvariable=self.list_ref,selectmode='extended', width = 30)
        button_select_sam = ttk.Button(self.window1,text='use selection as sample',command=self.select_sam)
        button_select_ref = ttk.Button(self.window1,text='use selection as reference',command=self.select_ref)
        button_calculate = ttk.Button(self.window1,text='calculate',command=self.calculate)
        button_select_combine_ref = ttk.Button(self.window1,text='select and combine',command=self.combine_select_ref)
        button_checktrans = ttk.Button(self.window1,text='check result',command=self.checktrans)
        button_trans_export = ttk.Button(self.window1, text = 'export transmission', command = self.export_complex_transmission)
        label_status = ttk.Label(self.window1,text='status: ')
        label_action = ttk.Label(self.window1,textvariable=self.action)
        label_select_sam = ttk.Label(self.window1, text = 'Selected Sample Specs: \n')
        label_select_ref = ttk.Label(self.window1, text = 'Selected Reference Specs: \n')
        label_show_selected_sam = ttk.Label(self.window1, textvariable = self.sam_selected)
        label_show_selected_ref = ttk.Label(self.window1, textvariable = self.ref_selected)
        
        
        self.listbox_sam.grid(column=0,row=0, columnspan = 2, sticky = 'w')
        self.listbox_ref.grid(column=2,row=0,columnspan = 2, sticky = 'w') 
        label_select_sam.grid(column = 0, row = 1, sticky = 'w')
        label_show_selected_sam.grid(column = 1, row = 1, sticky = 'w')
        label_select_ref.grid(column = 2, row = 1, sticky = 'w')
        label_show_selected_ref.grid(column = 3, row = 1, sticky = 'w')
        
        button_select_sam.grid(column=0,row=2, sticky = 'w')
        button_select_ref.grid(column=2,row=2, sticky = 'w')
        button_select_combine_ref.grid(column=3,row=2, sticky = 'w')
        button_calculate.grid(column=4,row=3, sticky = 'w')
        button_checktrans.grid(column=0,row=4, sticky = 'w')
        button_trans_export.grid(column = 1, row = 4, sticky = 'w')
        label_status.grid(column=0,row=5, sticky = 'w')
        label_action.grid(column=1,row=5, sticky = 'w')
        
        
        
        # self.window2.grid(column=1,row=0,sticky='nsew')
        
        
    def select_sam(self):
        self.idx_sam = self.listbox_sam.curselection()
        display = ''
        for i in self.idx_sam:
            display = display + self.fnames_sam[i] + '\n'
        self.sam_selected.set(display)    
        self.action.set('sample scans selected')
        
    def select_ref(self):
        self.idx_ref = self.listbox_ref.curselection()
        display = ''
        for i in self.idx_ref:
            display = display + self.fnames_ref[i] + '\n'
        self.ref_selected.set(display) 
        self.action.set('reference scans selected')
    
    def calculate(self):
        self.transmission_cx = {}
        self.transmission = {}
        self.freq_plot = {}
        if len(self.idx_sam) == len(self.idx_ref) and len(self.idx_ref) != 1:
            for i in range(0,len(self.idx_sam)):
               self.transmission[self.fnames_sam[self.idx_sam[i]]] = abs(self.sxsam[self.fnames_sam[self.idx_sam[i]]])/abs(self.sxref[self.fnames_ref[self.idx_ref[i]]])
               self.transmission_cx[self.fnames_sam[self.idx_sam[i]]] = self.sxsam[self.fnames_sam[self.idx_sam[i]]]/self.sxref[self.fnames_ref[self.idx_ref[i]]]
               self.freq_plot[self.fnames_sam[self.idx_sam[i]]] = self.freq[self.fnames_sam[self.idx_sam[i]]]
        
        if len(self.idx_ref) == 1 and len(self.idx_sam) == 1:
            self.transmission[self.fnames_sam[self.idx_sam[0]]] = abs(self.sxsam[self.fnames_sam[self.idx_sam[0]]])/abs(self.sxref[self.fnames_ref[self.idx_ref[0]]])
            self.transmission_cx[self.fnames_sam[self.idx_sam[0]]] = self.sxsam[self.fnames_sam[self.idx_sam[0]]]/self.sxref[self.fnames_ref[self.idx_ref[0]]]
            self.freq_plot[self.fnames_sam[self.idx_sam[0]]] = self.freq[self.fnames_sam[self.idx_sam[0]]]
        self.action.set('transmission calculated')
        
    def export_complex_transmission(self):
        path = tk.filedialog.askdirectory()
        for key in list(self.transmission_cx):
            filename = path + '/' + 'trans_' + key + '_1.dat'
            # trans_real = np.real(self.transmission_cx[key])
            # trans_imag = np.imag(self.transmission_cx[key])
            # trans_out = np.stack((trans_real, trans_imag), axis = 1)
            # np.savetxt(filename, trans_out)
            trans_out = np.stack((self.freq_plot[key], self.transmission_cx[key]), axis = 1)
            np.savetxt(filename, trans_out)
                
    def combine_select_ref(self):
        self.idx_ref = self.listbox_datas.curselection()
        for i in self.idx_ref:
            self.sxref[self.fnames_ref[i]] = (self.sxref[self.fnames_ref[i]] + self.sxref[self.fnames_ref[i+1]])/2
        self.action.set('reference combined')
            
    def checktrans(self):
        self.window2 = tk.Toplevel()
        self.st_plotgo = pt.Plottool(self.window2, self.freq_plot, self.transmission)


class study_index:
    def __init__(self,frame,x,t, path = None):
        self.x = copy.deepcopy(x)
        self.t = copy.deepcopy(t)
        self.x_backup = copy.deepcopy(x)
        self.t_backup = copy.deepcopy(t)
        if path != None:
            self.path = path
        # self.x = dp.Basic_functions().normalize_signal_samref(self.x)
        self.frame_si = frame
        self.window_1 = ttk.Labelframe(self.frame_si,text='main')
        self.window_2 = ttk.Labelframe(self.frame_si,text='display')
        self.window_3 = ttk.LabelFrame(self.frame_si,text='status')
        self.window_4 = ttk.Labelframe(self.frame_si, text = 'save results')
        self.window_1.grid(column=0,row=0)
        self.window_2.grid(column=1,row=0)
        self.window_3.grid(column=0,row=2)
        self.window_4.grid(column = 0, row = 1)
        self.fnames = list(x)
        self.list = tk.StringVar(value=self.fnames)
        # self.si_action = tk.StringVar(value='ready')
        self.si_listbox= tk.Listbox(self.window_1,listvariable=self.list,selectmode='extended', width = 50)
        self.si_idx_echo = None
        try:
            open('study.ini')
        except:
            self.si_action = tk.StringVar(value = 'Ready with no default values loaded') 
            self.choppoints_sam = tk.StringVar()
            self.choppoints_ref = tk.StringVar()
            self.choppoints_echo = tk.StringVar()
            self.L = tk.DoubleVar(value = 0.00108)
        else: 
            config = cfp.ConfigParser()
            config.read('study.ini')
            self.si_action = tk.StringVar(value = 'Ready with default values loaded') 
            self.choppoints_sam = tk.StringVar(value = config['si_default_values']['Choppoints_sam'])
            self.choppoints_ref = tk.StringVar(value = config['si_default_values']['Choppoints_ref'])
            self.choppoints_echo = tk.StringVar(value = config['si_default_values']['Choppoints_echo'])
            self.L = tk.DoubleVar(value = float(config['si_default_values']['Sample_thickness']))
            
        self.phasemod = tk.IntVar(value = 0)
        self.f0 = tk.DoubleVar(value = 0.25)
        self.f1 = tk.DoubleVar(value = 1.6)
        # self.index_save = {}
        # self.absorption_save = {}
        self.optical_constants = {}
        self.freq_save = {}
        
        
        '''window content for window_1: operate spectrums'''
        
        
        button_select_sam = ttk.Button(self.window_1,text='use selection as sample',command=self.si_select_sam)
        button_select_ref = ttk.Button(self.window_1,text='use selection as reference',command=self.si_select_ref)
        button_select_echo = ttk.Button(self.window_1,text='use selection as echo',command=self.si_select_echo)
        button_calculate = ttk.Button(self.window_1,text='calculate',command=self.calculate_index)
        button_show = ttk.Button(self.window_1,text='display result',command=self.display)
        
        label_L = ttk.Label(self.window_1,text='Enter estimated thickness: ')
        entry_L = ttk.Entry(self.window_1,textvariable=self.L,width=10)
        label_choppoints_sam = ttk.Label(self.window_1, text='Enter times to chop sample (start, stop): ')
        entry_choppoints_sam = ttk.Entry(self.window_1, textvariable = self.choppoints_sam, width=10)
        label_choppoints_ref = ttk.Label(self.window_1, text='Enter times to chop ref (start, stop): ')
        entry_choppoints_ref = ttk.Entry(self.window_1, textvariable = self.choppoints_ref, width=10)
        label_choppoints_echo = ttk.Label(self.window_1, text='Enter times to build echo from sample (start, stop): ')
        entry_choppoints_echo = ttk.Entry(self.window_1, textvariable = self.choppoints_echo, width=10)
        button_chopspec_sam = ttk.Button(self.window_1,text='chop sample spec',command=self.chop_sam_spec)
        button_chopspec_ref = ttk.Button(self.window_1,text='chop reference spec',command=self.chop_ref_spec)
        button_buildpec_echo = ttk.Button(self.window_1,text='chop echo spec',command=self.build_echo_spec)
        button_padspec = ttk.Button(self.window_1,text='pad specs', command=self.pad_spec)
        button_save_default = ttk.Button(self.window_1, text = 'save as default', command = self.build_config_file)
        button_show_selected = ttk.Button(self.window_1, text = 'show selected', command = self.plot_selected)
        button_restore = ttk.Button(self.window_1, text = 'restore spectrums', command = self.restore_backups)
        
        checkbutton_phasemod = ttk.Checkbutton(self.window_1, text = 'add 2PI to sam phase', variable = self.phasemod)
        
        self.si_listbox.grid(column=0,row=0,columnspan=7, sticky = 'w')
        button_select_sam.grid(column=0,row=1)
        button_select_ref.grid(column=1,row=1)
        button_select_echo.grid(column=2,row=1)
        button_show_selected.grid(column = 3, row = 1)
        label_L.grid(column=0,row=2)
        entry_L.grid(column=1,row=2)
        button_calculate.grid(column=2,row=3)
        button_show.grid(column=3,row=3)
        label_choppoints_sam.grid(column = 0,row = 4)
        entry_choppoints_sam.grid(column = 1,row = 4)
        button_chopspec_sam.grid(column = 2, row = 4)
        button_padspec.grid(column = 3,row = 4)
        checkbutton_phasemod.grid(column = 4, row = 4, sticky = 'w')
        label_choppoints_ref.grid(column = 0,row = 5)
        entry_choppoints_ref.grid(column = 1,row = 5)
        button_chopspec_ref.grid(column = 2, row = 5)
        label_choppoints_echo.grid(column = 0,row = 6)
        entry_choppoints_echo.grid(column = 1,row = 6)
        button_buildpec_echo.grid(column = 2, row = 6)
        button_save_default.grid(column = 0, row = 7)
        button_restore.grid(column = 0, row = 8)
        
        
        
        
        '''window content for save results window'''
        
        # button_save = ttk.Button(self.window_4, text = 'save results', command = self.save)
        button_plotsave = ttk.Button(self.window_4, text = 'plot saved results', command = self.Edit_plot_save)
        # entry_f0 = ttk.Entry(self.window_4, textvariable = self.f0, width = 10)
        # entry_f1 = ttk.Entry(self.window_4, textvariable = self.f1, width = 10)
        # label_f0 = ttk.Label(self.window_4, text = 'Freq start from: ')
        # label_f1 = ttk.Label(self.window_4, text = 'Freq stop at: ')
        # label_f0.grid(column = 0, row = 0)
        # entry_f0.grid(column = 1, row = 0)
        # label_f1.grid(column = 2, row = 0)
        # entry_f1.grid(column = 3, row = 0)
        # button_save.grid(column = 0, row = 1)
        button_plotsave.grid(column = 0, row = 0)
        # button_save.grid(column = 0, row = 5)
        self.fig, self.ax = plt.subplots(1,1, figsize = (16,9))
        self.ax2 = self.ax.twinx()
        self.canvas0 = FigureCanvasTkAgg(figure = self.fig, master = self.window_2)
        self.canvas0.get_tk_widget().grid(column=0,row=0,columnspan=6)
        button_examplot = ttk.Button(self.window_2,text='exam plot',command=self.exam_plot)
        button_examplot.grid(column=0,row=1)
        
        '''window for status bar'''
        
        
        label_action = ttk.Label(self.window_3, textvariable = self.si_action)
        label_action.grid(column=0,row=0)
        
        self.bf = dp.Basic_functions()
    
    def si_select_sam(self):
        self.si_idx_sam = self.si_listbox.curselection()[0]
        self.si_action.set('sample spectrum selected')
        
    def si_select_ref(self):
        self.si_idx_ref = self.si_listbox.curselection()[0]
        self.si_action.set('reference spectrum selected')
        
    def si_select_echo(self):
        self.si_idx_echo = self.si_listbox.curselection()[0]
        self.si_action.set('echo spectrum selected')
    
    def restore_backups(self):
        self.x = copy.deepcopy(self.x_backup)
        self.t = copy.deepcopy(self.t_backup)
    
    def build_config_file(self):
        si_config = cfp.ConfigParser()
        si_config['si_default_values'] = {'Sample_thickness': str(self.L.get()),
                                           'Choppoints_sam': self.choppoints_sam.get(),
                                           'Choppoints_ref': self.choppoints_ref.get(),
                                           'Choppoints_echo': self.choppoints_echo.get()}
        with open('study.ini', 'w') as si_default_values:
            si_config.write(si_default_values)
        self.si_action.set('Current entry values are saved as default values')
            
    def plot_selected(self):
        fig, ax = plt.subplots(1, 1, figsize = (16,9))
        self.canvas0.figure = fig
        key_sam = self.fnames[self.si_idx_sam]
        key_ref = self.fnames[self.si_idx_ref]
        ax.plot(self.t[key_sam], self.x[key_sam], label = 'sample')
        ax.plot(self.t[key_ref], self.x[key_ref], label = 'reference')
        ax.legend()
        self.canvas0.draw()
        
    def calculate_index(self):  
        # f_bot = self.f0.get()
        # self.x = dp.Basic_functions().weak_signal_remove(self.x)
        freq,sx = dp.Basic_functions().fftx_hilbert(self.x, self.t, 2)
        freq_sam = freq[self.fnames[self.si_idx_sam]]
        # freq_low = freq_sam[np.where(freq_sam <= 2)]
        # self.freq_sam = freq_sam
        sx_sam = sx[self.fnames[self.si_idx_sam]]
        # sx_sam = sx_sam[np.where(freq_sam <= 2)]
        # freq_ref = freq[self.fnames[self.si_idx_ref]]
        sx_ref = sx[self.fnames[self.si_idx_ref]]
        # sx_ref = sx_ref[np.where(freq_sam <= 2)]
        if self.phasemod.get() == 0:
            if self.si_idx_echo == None:
                freq_new, sx_sam_new, sx_ref_new, self.index, self.absorption, self.phi_sam, self.phi_ref = self.bf.getindex_array(freq_sam, sx_sam,
                                                                                                                 sx_ref,self.L.get())
            else:
                # freq_echo = freq[self.fnames[self.si_idx_echo]]
                sx_echo = sx[self.fnames[self.si_idx_echo]]
                sx_echo = sx_echo[np.where(freq_sam <= 2)]
                freq_new, sx_sam_new, sx_ref_new, self.index, self.absorption, self.phi_sam, self.phi_ref = self.bf.getindex_array(freq_sam, sx_sam, 
                                                                                                                 sx_ref, self.L.get(),
                                                                                                                 SX_echo=sx_echo)
        else:
            if self.si_idx_echo == None:
                freq_new, sx_sam_new, sx_ref_new, self.index, self.absorption, self.phi_sam, self.phi_ref = self.bf.getindex_array(freq_sam, sx_sam,
                                                                                                                 sx_ref,self.L.get(), ph_mod = 2*pi)
            else:
                # freq_echo = freq[self.fnames[self.si_idx_echo]]
                sx_echo = sx[self.fnames[self.si_idx_echo]]
                sx_echo = sx_echo[np.where(freq_sam <= 2)]
                freq_new, sx_sam_new, sx_ref_new, self.index, self.absorption, self.phi_sam, self.phi_ref = self.bf.getindex_array(freq_sam, sx_sam, 
                                                                                                                 sx_ref, self.L.get(),
                                                                                                                 SX_echo=sx_echo, ph_mod = 2*pi)
        self.freq = freq_new
        # self.canvas0.figure = Figure(figsize=(15,8))
        # self.ax.cla()
        # plot0 = self.ax
      
        # plot0.plot(freq_new, abs(sx_sam_new), freq_new, abs(sx_ref_new))
        
        # if self.si_idx_echo != None:
        #     plot0.plot(freq_new,abs(sx_echo))
      
        # plot0.set_xlabel('frequency (THz)',fontsize=15)
        # plot0.set_ylabel('amplitude (V)',fontsize=15)
        # plot0.grid(1)
        # self.canvas0.draw()
        
        self.freq_save[self.fnames[self.si_idx_sam]] = self.freq
        self.optical_constants[self.fnames[self.si_idx_sam]+'_index'] = self.index
        self.optical_constants[self.fnames[self.si_idx_sam]+'_absorb'] = self.absorption
        
        dp.Basic_functions.save_data(self.path+self.fnames[self.si_idx_sam]+'_index_1.dat', self.freq, self.index)
        dp.Basic_functions.save_data(self.path+self.fnames[self.si_idx_sam]+'_absorption_1.dat', self.freq, self.absorption)
        # self.index_save[self.fnames[self.si_idx_sam]] = self.index
        # self.absorption_save[self.fnames[self.si_idx_sam]] = self.absorption
        
        self.si_action.set('computation complete')
    
    
    def display(self):
        self.canvas0.figure,self.ax = plt.subplots(1,1,figsize=(15,8)) 
        self.ax2 = self.ax.twinx()
        self.ax.cla()
        self.ax2.cla()
        self.ax.plot(self.freq,self.index, color = 'red', label = self.fnames[self.si_idx_sam])
        self.ax.set_xlim(left=0,right=1.65)
        self.ax.set_xlabel('frequency (THz)',fontsize=15)
        self.ax.set_ylabel('refractive index',fontsize=15, color  = 'red')
        self.ax.legend(loc = 'best', fontsize = 15)
        # self.ax.set_ylim(bottom = 0, top = 2)a.
        self.ax2.plot(self.freq,self.absorption, color = 'blue')
        self.ax2.set_xlim(left=0,right=1.65)
        self.ax2.set_xlabel('frequency (THz)',fontsize=15)
        self.ax2.set_ylabel('absorption',fontsize=15, color = 'blue')
        self.canvas0.draw()
        
        
    def save(self):
        self.freq_save[self.fnames[self.si_idx_sam]] = self.freq
        self.optical_constants[self.fnames[self.si_idx_sam]+'_index'] = self.index
        self.optical_constants[self.fnames[self.si_idx_sam]+'_absorb'] = self.absorption
        # self.index_save[self.fnames[self.si_idx_sam]] = self.index
        # self.absorption_save[self.fnames[self.si_idx_sam]] = self.absorption
        self.si_action.set('result saved')
        
    # def plot_save(self):
    #     self.canvas0.figure, self.ax = plt.subplots(1,1,figsize=(15,8)) 

    #     for key in list(self.index_save):
    #         self.ax[0].plot(self.freq_save[key],self.index_save[key], label = key)
    #         self.ax[1].plot(self.freq_save[key],self.absorption_save[key], label = key)
    #     # self.ax[0].set_xlim(left=0,right=2)
    #     self.ax[0].set_xlabel('frequency (THz)',fontsize=15)
    #     self.ax[0].set_ylabel('refractive index',fontsize=15)
    #     self.ax[0].legend(loc = 'best', fontsize = 15)
    #     # self.ax[0].set_ylim(bottom = 0, top = 2)
    #     # self.ax[1].set_xlim(left=0,right=2)
    #     self.ax[1].set_xlabel('frequency (THz)',fontsize=15)
    #     self.ax[1].set_ylabel('absorption',fontsize=15)
    #     self.ax[1].legend(loc = 'best', fontsize = 15)
    #     self.canvas0.draw()
    #     self.si_action.set('Figure plotted')
        
    def Edit_plot_save(self):
        # self.canvas0.figure,self.ax = plt.subplots(2,1,figsize=(15,8)) 
        window_temp1 = tk.Toplevel()
        plotgo = pt.Plottool(window_temp1, self.freq_save, self.optical_constants, nrow = 1, ncol = 1)
        
        
    def exam_plot(self):
        root = tk.Tk()
        frame = ttk.Frame(root)
        frame.pack(fill=tk.BOTH,expand=True)
        self.canvas1 = FigureCanvasTkAgg(self.canvas0.figure,master=frame)
        self.canvas1.figure = self.canvas0.figure
        self.canvas1.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.canvas1.draw()
        self.toolbar = NavigationToolbar2Tk(self.canvas1,frame)
        
        
        
        
    def build_echo_spec(self):
        key_sam = self.fnames[self.si_idx_sam]
        t_points = self.choppoints_echo.get().split(',')
        p0 = float(t_points[0])
        p1 = float(t_points[1])
        idx0 = np.where(self.t[key_sam] >= p0)[0][0]
        idx1 = np.where(self.t[key_sam] <= p1)[0][-1]
        key_echo = key_sam + '_echo'
        t0 = self.t[key_sam][0 : idx0]
        # t1 = self.t[key][idx0:idx1]
        t1 = self.t[key_sam][idx1 : ]
        pad0 = np.zeros_like(t0)
        pad1 = np.zeros_like(t1)
        
        self.x[key_echo] = np.concatenate((pad0, self.x[key_sam][idx0:idx1], pad1), axis = 0)
        self.t[key_echo] = self.t[key_sam]
        self.fnames = list(self.x)
        self.list.set(self.fnames)
        self.si_action.set('Echo spectrum built')
        
        
    def chop_ref_spec(self):
        key = self.fnames[self.si_idx_ref]
        t_points = self.choppoints_ref.get().split(',')
        p0 = round(float(t_points[0]), ndigits = 2)
        p1 = round(float(t_points[1]), ndigits = 2)
        idx0 = np.where(self.t[key] >= p0)[0][0]
        idx1 = np.where(self.t[key] <= p1)[0][-1]
        t0 = self.t[key][0 : idx0]
        # t1 = self.t[key][idx0:idx1]
        t1 = self.t[key][idx1 : ]
        pad0 = np.zeros_like(t0)
        pad1 = np.zeros_like(t1)
        
        self.x[key] = np.concatenate((pad0, self.x[key][idx0:idx1], pad1), axis = 0)
        # self.t[key] = self.t[key][idx0:idx1]
        self.si_action.set('reference spectrum chopped')
        
    def chop_sam_spec(self):
        key = self.fnames[self.si_idx_sam]
        t_points = self.choppoints_sam.get().split(',')
        
        p0 = round(float(t_points[0]), ndigits = 2)
        p1 = round(float(t_points[1]), ndigits = 2)
        
        idx0 = np.where(self.t[key] >= p0)[0][0]
        idx1 = np.where(self.t[key] <= p1)[0][-1]
        
        t0 = self.t[key][0 : idx0]
        # t1 = self.t[key][idx0:idx1]
        t1 = self.t[key][idx1 : ]
        
        pad0 = np.zeros_like(t0)
        pad1 = np.zeros_like(t1)
        
        self.x[key] = np.concatenate((pad0, self.x[key][idx0:idx1], pad1), axis = 0)
        # self.t[key] = self.t[key][idx0:idx1]
        self.si_action.set('sample spectrum chopped')
        
    def pad_spec(self):
        key_sam = self.fnames[self.si_idx_sam]
        key_ref = self.fnames[self.si_idx_ref]
        t_sam = np.round(self.t[key_sam], 2)
        t_ref = np.round(self.t[key_ref], 2)
        # self.t_ref = t_ref
        # x_ref = self.t[key_ref]
        t_max = max(t_sam[-1], t_ref[-1])
        # self.t_max = t_max
        t_min = min(t_sam[0], t_ref[0])
        if self.si_idx_echo != None:
            key_echo = self.fnames[self.si_idx_echo]
            t_echo = np.round(self.t[key_echo], 2)
            t_max = max(t_max, t_echo[-1])
            t_min = min(t_min, t_echo[0])
        
        timestep = t_sam[1] - t_sam[0]
        t_all = np.round(np.arange(t_min - 5, t_max + 5, step = timestep), 2)
        self.t[key_sam] = t_all
        self.t[key_ref] = t_all
        
        self.t_sam_0 = t_all[np.where(t_all < t_sam[0])]
        self.t_sam_1 = t_all[np.where(t_all > t_sam[-1])]
        t_sam_0 = t_all[np.where(t_all < t_sam[0])]
        t_sam_1 = t_all[np.where(t_all > t_sam[-1])]
        pad_sam_0 = np.zeros_like(t_sam_0) + self.x[key_sam][0]
        pad_sam_1 = np.zeros_like(t_sam_1) + self.x[key_sam][-1]
        self.x[key_sam] = np.concatenate((pad_sam_0, self.x[key_sam], pad_sam_1), axis = 0)
        
        t_ref_0 = t_all[np.where(t_all < t_ref[0])]
        t_ref_1 = t_all[np.where(t_all > t_ref[-1])]
        pad_ref_0 = np.zeros_like(t_ref_0) + self.x[key_ref][0]
        pad_ref_1 = np.zeros_like(t_ref_1) + self.x[key_ref][-1]
        self.x[key_ref] = np.concatenate((pad_ref_0, self.x[key_ref], pad_ref_1), axis = 0)
        
        if self.si_idx_echo != None:
            t_echo_0 = t_all[np.where(t_all < t_echo[0])]
            t_echo_1 = t_all[np.where(t_all > t_echo[-1])]
            pad_echo_0 = np.zeros_like(t_echo_0) + self.x[key_echo][0]
            pad_echo_1 = np.zeros_like(t_echo_1) + self.x[key_echo][-1]
            self.x[key_echo] = np.concatenate((pad_echo_0, self.x[key_echo], pad_echo_1), axis = 0)
            self.t[key_echo] = t_all
        
        
            
        self.figure, self.ax = plt.subplots(3, 1, figsize = (16,9))
        self.canvas0.figure = self.figure
        
        plot1 = self.ax[0]
        plot2 = self.ax[1]
        plot3 = self.ax[2]
        plot1.plot(self.t[key_sam] , self.x[key_sam],label = 'sample')
        plot2.plot(self.t[key_ref] , self.x[key_ref],label = 'reference')
        if self.si_idx_echo != None:
            plot3.plot(self.t[key_echo], self.x[key_echo],label = 'echo')
            plot3.set_xlabel('time (ps)')
            plot3.set_ylabel('amplitude (V)')
            plot3.legend(loc='best')
        plot1.set_xlabel('time (ps)')
        plot1.set_ylabel('amplitude (V)')
        plot1.legend(loc='best')
        plot2.set_xlabel('time (ps)')
        plot2.set_ylabel('amplitude (V)')
        plot2.legend(loc='best')
       
        self.canvas0.draw()
        
        
        
class study_polarimetry:
    def __init__(self,frame,t, x, datapath = None):
        self.path = datapath
        self.frame_po = frame
        self.t = copy.deepcopy(t)
        self.x = copy.deepcopy(x)
        self.freq_x, self.sx = dp.Basic_functions().fftx(self.x, self.t, 2)
        self.sx_abs = dp.Basic_functions().dict_getabs(self.sx)

        self.polarangle = dict()
        self.ellipticity = dict()
        self.polarangle_reduced = {}
        self.ellipticity_reduced = {}
        self.transmission = {}
        self.transmission_normal = {}
        self.freq_plot = {}
        self.sort_keyword = tk.StringVar(value = 'T')
        
        '''first window for polarimetry calculation'''
        
        #variables
        self.fnames = list(x)
        self.list_name = tk.StringVar(value=self.fnames)
        self.selection_sam = tk.StringVar(value='Selected sample spectrums: \n')
        self.selection_ref = tk.StringVar(value='Selected reference spectrums: \n')
        self.action = tk.StringVar(value='ready')
        
        
        # mainwindow
        self.window1 = ttk.Labelframe(self.frame_po,text='main')
        self.window1.grid(column=0,row=0)
        self.choppoints_sam = tk.StringVar(value = '0,0')
        self.choppoints_ref = tk.StringVar(value = '0,0')
        
        
        # widges
        self.listbox= tk.Listbox(self.window1,listvariable = self.list_name,selectmode='extended', width = 50)
        button_sortlist = ttk.Button(self.window1,text='sort x and y lists', command = self.sort_input_list)
        cbox_sort_keyword = ttk.Combobox(self.window1, textvariable = self.sort_keyword, width = 10)
        cbox_sort_keyword['value'] = ['K', 'T']
        label_select_sam = ttk.Label(self.window1,textvariable = self.selection_sam)
        label_select_ref = ttk.Label(self.window1,textvariable = self.selection_ref)
        
        button_auto_select = ttk.Button(self.window1, text = 'auto select spec', command = self.autoselect)
        button_select_sam = ttk.Button(self.window1,text = 'use selection as +45',command=self.select_sam)
        button_select_ref = ttk.Button(self.window1,text = 'use selection as -45',command=self.select_ref)
        # button_calculate_transmission = ttk.Button(self.window1,text='calculate transmission',command=self.calculate_transmission)
        button_calculate_polarimetry = ttk.Button(self.window1,text='calculate polarimetry',
                                                  command = self.calculate_polarimetry)
        button_show_selected = ttk.Button(self.window1, text = 'show selected', command = self.plot_selected)
        button_checktrans = ttk.Button(self.window1,text='show transmission',command=self.checktrans)
        button_show_polarangle = ttk.Button(self.window1,text='show polarangle',command=self.show_polarangle)
        button_show_ellipticity = ttk.Button(self.window1,text='show ellipticity',command=self.show_ellipticity)
        button_show_polarangle_re = ttk.Button(self.window1,text='show reduced polarangle',command=self.show_polarangle_reduced)
        button_show_ellipticity_re = ttk.Button(self.window1,text='show reduced ellipticity',command=self.show_ellipticity_reduced)
        
        
        label_status = ttk.Label(self.window1,text='status: ')
        label_action = ttk.Label(self.window1,textvariable=self.action)
        
        
        self.listbox.grid(column=0,row=0,columnspan=2)
        label_select_sam.grid(column=2,row=0,columnspan=2, sticky = 'nw')
        label_select_ref.grid(column=4,row=0,columnspan=2, sticky = 'nw')
        
        cbox_sort_keyword.grid(column = 0, row = 1, sticky = 'w')
        button_sortlist.grid(column = 1, row = 1, sticky = 'w')
        button_select_sam.grid(column = 2,row = 1, sticky = 'w')
        button_select_ref.grid(column = 3,row = 1, sticky = 'w')
        button_auto_select.grid(column = 4, row = 1,sticky = 'w')
        button_show_selected.grid(column = 5, row = 1, sticky = 'w')
        button_calculate_polarimetry.grid(column=0,row=2)
        # button_calculate_transmission.grid(column=1,row=2)
        button_show_polarangle.grid(column = 0, row = 3, sticky = 'w')
        button_checktrans.grid(column = 1 ,row = 3, sticky = 'w')
        button_show_ellipticity.grid(column = 2, row = 3, sticky = 'w')
        button_show_polarangle_re.grid(column = 0, row = 4)
        button_show_ellipticity_re.grid(column = 1, row = 4)
        
        
        
        label_status.grid(column = 0, row = 5, sticky = 'w')
        label_action.grid(column = 1, row = 5, sticky = 'w')
        
        
        # self.window2 = ttk.Labelframe(self.frame_po, text='Edit Spectrums')
        # self.window2.grid(column = 0, row = 1, sticky='w')
        
        # label_choppoints_sam = ttk.Label(self.window2, text='Enter time points to chop sample (start, stop): ')
        # entry_choppoints_sam = ttk.Entry(self.window2, textvariable = self.choppoints_sam, width=10)
        # label_choppoints_ref = ttk.Label(self.window2, text='Enter time points to chop ref (start, stop): ')
        # entry_choppoints_ref = ttk.Entry(self.window2, textvariable = self.choppoints_ref, width=10)
        # button_chop = ttk.Button(self.window2, text = 'chop specs', command = self.spec_chop)
        
        # label_choppoints_sam.grid(column = 0, row = 0, sticky = 'w')
        # entry_choppoints_sam.grid(column = 1, row = 0, sticky = 'w')
        # label_choppoints_ref.grid(column = 2, row = 0, sticky = 'w')
        # entry_choppoints_ref.grid(column = 3, row = 0, sticky = 'w')
        # button_chop.grid(column = 0, row = 1, sticky = 'w')

        # self.window3 = ttk.LabelFrame(self.frame_po,text='plot')
        # self.window3.grid(column = 1,row = 0,sticky='w')
        # self.fig, self.ax = plt.subplots(2,2,figsize=(16,9))
        # self.canvas0 = FigureCanvasTkAgg(figure = self.fig, master = self.window3)
        # self.canvas0.get_tk_widget().grid(column=0,row=0,columnspan=6,rowspan=6,sticky='w')
        # self.canvas0.draw()
        
        # button_clear_plot = ttk.Button(self.window3, text = 'clear', command = self.plot_clear)
        
        # button_clear_plot.grid(column = 0, row = 7)
        
        self.window4 = ttk.Labelframe(self.frame_po, text = 'status')
        self.window4.grid(column = 0, row = 2)
        
        
    def autoselect(self):
        idx_sam = []
        idx_ref = []
        for i, name in enumerate(self.fnames):
            if 'p135' in name:
                idx_sam.append(i)
            if 'p45' in name:
                idx_ref.append(i)
        self.idx_sam = idx_sam
        self.idx_ref = idx_ref
        for i in self.idx_sam:
            self.selection_sam.set(self.selection_sam.get()+self.fnames[i]+'\n')
        for i in self.idx_ref:
            self.selection_ref.set(self.selection_ref.get()+self.fnames[i]+'\n')
        self.action.set('+45 and -45 spectrums selected')
        
                
        
        
    def select_sam(self):
        self.idx_sam = self.listbox.curselection()
        self.selection_sam.set('+45 deg spectrums: \n')
        for i in self.idx_sam:
            self.selection_sam.set(self.selection_sam.get()+self.fnames[i]+'\n')
        
        # for i in self.idx_sam:
        #     key = self.fnames[i]
        #     self.freq_plot[key] = self.freq_x[key]
        self.action.set('+45 deg scans selected')
        
    def select_ref(self):
        self.idx_ref = self.listbox.curselection()
        self.selection_ref.set('-45 deg spectrums: \n')
        for i in self.idx_ref:
            self.selection_ref.set(self.selection_ref.get()+self.fnames[i]+'\n')
        self.action.set('-45 deg scans selected')
    
    def plot_selected(self):
        for axes1, axes2 in self.ax:
            axes1.cla()
            axes2.cla()
        for i in self.idx_sam:
            key_sam = self.fnames[i]
            self.ax[0,0].plot(self.t[key_sam], self.x[key_sam], label = key_sam + '_x')
            self.ax[0,1].plot(self.t[key_sam], self.y[key_sam], label = key_sam + '_y')
        for j in self.idx_ref:
            key_ref = self.fnames[j]
            self.ax[1,0].plot(self.t[key_ref], self.x[key_ref], label = key_ref + '_x')
            self.ax[1,1].plot(self.t[key_ref], self.y[key_ref], label = key_ref + '_y')
        self.ax[0,0].legend(loc= 'best')
        self.ax[0,1].legend(loc= 'best')
        self.ax[1,0].legend(loc= 'best')
        self.ax[1,1].legend(loc= 'best')
        
        self.canvas0.draw()
        
    def plot_clear(self):
        for axes1, axes2 in self.ax:
            axes1.cla()
            axes2.cla()
        self.canvas0.draw()
        
    def sort_input_list(self):
        if self.sort_keyword.get() == 'K':
            self.fnames.sort(key = bf.dict_key_get_T)
            self.list_name.set(self.fnames)
        if self.sort_keyword.get() == 'T':
            self.fnames.sort(key = bf.dict_key_get_B)
            self.list_name.set(self.fnames)
            
    # def calculate_transmission(self):
    #     self.transmission_x = dict()
    #     self.transmission_y = dict()
    #     if len(self.idx_sam) == len(self.idx_ref) and len(self.idx_ref) != 1:
    #         for i in self.idx_sam:
    #             for j in self.idx_ref:
    #                 if self.fnames[i] in list(self.sy):
    #                     self.transmission_x[self.fnames[i]] = abs(self.sx[self.fnames[i]])/abs(self.sx[self.fnames[j]])
    #                     self.transmission_y[self.fnames[i]] = abs(self.sy[self.fnames[i]])/abs(self.sy[self.fnames[j]])
    #     if len(self.idx_ref) == 1:
    #         for i in self.idx_sam:
    #             self.transmission_x[self.fnames[i]] = abs(self.sx[self.fnames[i]])/abs(self.sx[self.fnames[self.idx_ref[0]]])
    #             self.transmission_y[self.fnames[i]] = abs(self.sy[self.fnames[i]])/abs(self.sy[self.fnames[self.idx_ref[0]]])
    #     self.action.set('transmission calculated')
     
    def spec_chop(self):
        for i in self.idx_sam:
            key_sam = self.fnames[i]
            t_points_sam = self.choppoints_sam.get().split(',')
            p0_sam = round(float(t_points_sam[0]), ndigits = 2)
            p1_sam = round(float(t_points_sam[1]), ndigits = 2)
            self.p0_sam = p0_sam
            self.p1_sam = p1_sam
            if p0_sam >= p1_sam:
                self.action.set('Error! Enter valid start and stop times for both sample and reference spectrums')
            else:
                self.x[key_sam] = dp.Basic_functions().array_chop_pad(self.t[key_sam], self.x[key_sam], p0_sam, p1_sam)
                self.y[key_sam] = dp.Basic_functions().array_chop_pad(self.t[key_sam], self.y[key_sam], p0_sam, p1_sam)
                
            # self.t[key] = self.t[key][idx0:idx1]
            
        for j in self.idx_ref:
            key_ref = self.fnames[j]
            t_points_ref = self.choppoints_ref.get().split(',')
            p0_ref = round(float(t_points_ref[0]), ndigits = 2)
            p1_ref = round(float(t_points_ref[1]), ndigits = 2)
            if p0_ref >= p1_ref:
                self.action.set('Error! Enter valid start and stop times for both sample and reference spectrums')
            else:
                self.x[key_ref] = dp.Basic_functions().array_chop_pad(self.t[key_ref], self.x[key_ref], p0_ref, p1_ref)
                self.y[key_ref] = dp.Basic_functions().array_chop_pad(self.t[key_ref], self.y[key_ref], p0_ref, p1_ref)
                self.action.set('sample and reference spectrum chopped')
        
    # def calculate_transmission(self):
    #     self.transmission_x = dict()
    #     self.transmission_y = dict()
    #     self.transmission = {}
    #     if len(self.idx_sam) == len(self.idx_ref):
    #         for i, j in zip(self.idx_sam, self.idx_ref):
    #             self.transmission_x[self.fnames[i]] = abs(self.sx[self.fnames[i]])/abs(self.sx[self.fnames[j]])
    #             self.transmission_y[self.fnames[i]] = abs(self.sy[self.fnames[i]])/abs(self.sy[self.fnames[j]]) 
    #             self.transmission[self.fnames[i]] = np.sqrt(self.transmission_x[self.fnames[i]]**2 + self.transmission_y[self.fnames[j]]**2)
    #     else:
    #         if len(self.idx_sam) < len(self.idx_ref) and 2*len(self.idx_sam) > len(self.idx_ref):  
    #             for i, value in enumerate(self.idx_sam):
    #                 key_sam = self.fnames[value]
    #                 key_ref1 = self.fnames[self.idx_ref[i]]
    #                 key_ref2 = self.fnames[self.idx_ref[i+1]]
    #                 self.transmission_x[key_sam] = abs(self.sx[key_sam])/abs((self.sx[key_ref1]+self.sx[key_ref2])/2)
    #                 self.transmission_y[key_sam] = abs(self.sy[key_sam])/abs((self.sy[key_ref1]+self.sy[key_ref2])/2)
    #                 self.transmission[key_sam] = np.sqrt(self.transmission_x[key_sam]**2 + self.transmission_y[key_sam]**2)
    #     self.action.set('transmission calculated')
        
    def calculate_polarimetry(self):  
        for i,j in zip(self.idx_sam, self.idx_ref):
            Ex = (self.sx[self.fnames[i]]-self.sx[self.fnames[j]])/np.sqrt(2)
            Ey = (self.sx[self.fnames[i]]+self.sx[self.fnames[j]])/np.sqrt(2)
            key_save = dp.Basic_functions().find_str_common(self.fnames[i], self.fnames[j])
            self.freq_plot[key_save] = self.freq_x[self.fnames[i]]
            Ecra = Ex+1j*Ey
            Ecri = Ex-1j*Ey
            self.transmission[key_save] = np.sqrt(abs(Ex)**2 + abs(Ey)**2)
            self.ellipticity[key_save] = (abs(Ecri)-abs(Ecra))/(abs(Ecri)+abs(Ecra))
            self.polarangle[key_save] = ((np.unwrap(np.angle(Ecra))-np.unwrap(np.angle(Ecri)))/2)/2/pi*360
        self.polarangle_phase_adjust()
        
        key0 = list(self.transmission)[0]
        for i, key in enumerate(list(self.transmission)):
            if i == 0:
                trans_0T = self.transmission[key]
            else:
                self.transmission_normal[key] = self.transmission[key]/trans_0T
        for key1, key in zip(list(self.freq_plot), list(self.polarangle)):
            fname_save = self.path + 'ellip'+key1 + '_1.dat'
            self.polarangle_reduced[key] = self.polarangle[key] - self.polarangle[key0]
            self.ellipticity_reduced[key] = self.ellipticity[key] - self.ellipticity[key0]
            ellip_save = np.stack((self.freq_plot[key1], self.ellipticity_reduced[key]), axis = 1)
            np.savetxt(fname_save, ellip_save)
        self.action.set('polarimetry calculated')
    
        
    def polarangle_phase_adjust(self):
        key = list(self.freq_plot)[0]
        idx0 = np.where(self.freq_plot[key]>0.1)[0][0]
        idx1 = np.where(self.freq_plot[key]<1.2)[0][-1]
        for i, key in enumerate(list(self.polarangle)):
            if i == 0:
                mean0 = np.mean(self.polarangle[key][idx0:idx1])
                while mean0 <= 0:
                    self.polarangle[key] = self.polarangle[key] + 180
                    mean0 = np.mean(self.polarangle[key][idx0:idx1])
            else:
                mean1 = np.mean(self.polarangle[key][idx0:idx1])
                # self.polarangle[key] = self.polarangle[key] - (mean1-mean0)
                
                while mean1 - mean0 >= 150:
                    self.polarangle[key] = self.polarangle[key] - 180
                    mean1 = np.mean(self.polarangle[key][idx0:idx1])
                while mean1 - mean0 <= -150:
                    self.polarangle[key] = self.polarangle[key] + 180
                    mean1 = np.mean(self.polarangle[key][idx0:idx1])
    
            
        
    def show_polarangle_reduced(self):
        newwindow = tk.Toplevel()
        pt.Plottool(newwindow, self.freq_plot, self.polarangle_reduced)
    
    def combine_select_ref(self):
        self.idx_ref = self.listbox_datas.curselection()
        for i in self.idx_ref:
            self.sxref[self.fnames_ref[i]] = (self.sxref[self.fnames_ref[i]] + self.sxref[self.fnames_ref[i+1]])/2
        self.action.set('reference combined')
            
    def checktrans(self):
        self.new_window = tk.Toplevel()
        pt.Plottool(self.new_window, self.freq_plot, self.transmission_normal)
        
    def show_polarangle(self):
        newwindow = tk.Toplevel()
        pt.Plottool(newwindow, self.freq_plot, self.polarangle)
    def show_ellipticity_reduced(self):
        newwindow = tk.Toplevel()
        pt.Plottool(newwindow, self.freq_plot, self.ellipticity_reduced)
    def show_ellipticity(self):
        newwindow = tk.Toplevel()
        pt.Plottool(newwindow, self.freq_plot, self.ellipticity)
        
        
# if __name__ == '__main__':
#     t = np
#     y = x**2
#     root=tk.Tk()
#     mainframe = ttk.Frame(root)
#     mainframe.grid(column = 0, row = 0)
#     pt=Plottool(mainframe, x, y)
#     root.mainloop()