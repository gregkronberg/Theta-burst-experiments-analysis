# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 11:57:16 2017

@author: Greg Kronberg

Preprocessing raw data files for LTP + DCS experiments
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.io as sio
import os

# directories for processed and raw files
fpath_r = 'D:\\Google Drive\\Work\\Research Projects\\Theta LTP\\Raw Matlab Data\\'
fpath_p = "D:\\Google Drive\\Work\\Research Projects\\Theta LTP\\Processed Matlab Data\\"

# list files in each directory
dir_r = os.listdir(fpath_r)
dir_p = os.listdir(fpath_p)

# list new unprocessed data files
dir_new = list(set(dir_r)-set(dir_p))

# experimental conditions
rec_loc = ['apical','basal','perforant','soma']
stim_loc = ['apical','basal','perforant']
stim_dcs = ['control','cathodal','anodal']
drug = ['none']

# organization of comments
com_columns= ['channel','block','sample','unknown','comment number']

# loop over new files
for slice in dir_new:
    matfile = sio.loadmat(fpath_r+slice)
    
    # extract experiment info from comments
    # loop over comments
    for idx,comment in enumerate(matfile['comtext']):
        com_chan = matfile['com'][idx,com_columns.index('channel')]
        
        # recording location
        if 'Soma' in comment:
            chan_soma = com_chan
        if 'Apical' in comment:
            chan_apical = com_chan
        if 'Basal' in comment:
            chan_basal = com_chan
        if 'Perforant' in comment:
            chan_perforant = com_chan
    
    
    