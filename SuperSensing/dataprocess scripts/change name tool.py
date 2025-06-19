# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:35:57 2023

@author: Admin
"""

import os
import tkinter as tk
from tkinter import ttk
import glob

class Simple_rename_tool:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("file name change tool")
        window = tk.LabelFrame(self.root)
        window.grid(column = 0, row = 0)
        self.oldname = tk.StringVar()
        self.newname = tk.StringVar()
        
        button_directory = ttk.Button(window, text = 'select folder', command = self.get_file_path)
        label_oldname = ttk.Label(window, text = 'Enter the name to be changed: ')
        entry_oldname = ttk.Entry(window, textvariable = self.oldname)
        label_newname = ttk.Label(window, text = 'Enter the name to change to: ')
        entry_newname = ttk.Entry(window, textvariable = self.newname)
        button_rename = ttk.Button(window, text = 'rename', command = self.change_file_names)
        
        button_directory.grid(column = 0, row = 0, sticky = 'w')
        label_oldname.grid(column = 0, row = 1, sticky = 'w')
        entry_oldname.grid(column = 1, row = 1, sticky = 'w')
        label_newname.grid(column = 0, row = 2, sticky = 'w')
        entry_newname.grid(column = 1, row = 2, sticky = 'w')
        button_rename.grid(column = 0, row = 3, sticky = 'w')
        
    def get_file_path(self):
        self.path = tk.filedialog.askdirectory()
        
    
    def change_file_names(self):
        all_fnames = glob.glob(self.path)
        for fname in all_fnames:
            if self.oldname.get() in fname:
                idx = fname.index(self.oldname.get())
                fname_new = fname[0:idx]+self.newname.get()+fname[idx+len(self.oldname.get()):-1]
                os.rename(fname, fname_new)
    
    
if __name__ == '__main__':
    srt = Simple_rename_tool()
    srt.root.mainloop()