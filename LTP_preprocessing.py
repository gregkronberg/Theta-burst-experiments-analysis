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

for slice in dir_new:
    matfile = sio.loadmat(fpath_r+slice)
    
    