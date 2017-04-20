# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 16:06:46 2017

@author: Greg Kronberg

Analyze area under fEPSP's during induction as predictor of synaptic plasticity
"""
#==============================================================================
# set up shell and modules
#==============================================================================
# reset shell
from IPython import get_ipython
get_ipython().magic('reset -sf')

# import modules
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
import scipy.io as sio
import os
import re 
import pylab

# close figures
plt.close('all')

#==============================================================================
# Set up directories
#==============================================================================
# choose computer to run on (0 = laptop, 1 = lab desktop)
comp = 1
if comp == 0: # laptop
    fpath_r = 'C:\\Users\\Greg Kronberg\\Google Drive\\Work\\Research Projects\
\\Theta LTP\\Raw Matlab Data'
    fpath_p = 'C:\Users\\Greg Kronberg\\Google Drive\\Work\\Research Projects\
\\Theta LTP\\Processed Python Data'
elif comp==1: # desktop
    fpath_r = 'D:\\Google Drive\\Work\\Research Projects\\Theta LTP\
\\Raw Matlab Data\\'
    fpath_p = 'D:\\Google Drive\\Work\\Research Projects\\Theta LTP\
\\Processed Python Data\\'

# list files in each directory
dir_r = os.listdir(fpath_r)
dir_p = os.listdir(fpath_p)

#==============================================================================
# experimental conditions
#==============================================================================
rec_loc = ['Apical','Basal','Perforant','Soma']
stim_loc = ['apical','basal','perforant']
stim_dcs = ['control','cathodal','anodal']
ind_num = ['induction_','induction2_','induction3_','induction4_']
slice_info = ['age_','hemi_','height_','current_']
drug = ['none']

# baseline recording period (minutes)
baset_pre = 20
baset_post = 60

#%%==============================================================================
# loop over new files
#==============================================================================
for slice in dir_new:
    # load matlab file
    matfile = sio.loadmat(fpath_p+slice)
    
    
