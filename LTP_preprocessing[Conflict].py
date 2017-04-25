# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 11:57:16 2017

@author: Greg Kronberg

Preprocessing raw data files for LTP + DCS experiments
"""

import numpy as np
import scipy as sp
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.io as sio
import os


#==============================================================================
# # directories for processed and raw files
#==============================================================================
# laptop
fpath_r = 'C:\\Users\\Greg Kronberg\\Google Drive\Work\\Research Projects\\Theta LTP\\Raw Matlab Data\\'
fpath_p = 'C:\\Users\\Greg Kronberg\\Google Drive\Work\\Research Projects\\Theta LTP\\Processed Matlab Data\\'

# desktop
#fpath_r = 'D:\\Google Drive\\Work\\Research Projects\\Theta LTP\\Raw Matlab Data\\'
#fpath_p = 'D:\\Google Drive\\Work\\Research Projects\\Theta LTP\\Processed Matlab Data\\'

# list files in each directory
dir_r = os.listdir(fpath_r)
dir_p = os.listdir(fpath_p)

# list new unprocessed data files
dir_new = list(set(dir_r)-set(dir_p))

#==============================================================================
#  experimental conditions
#==============================================================================
rec_loc = ['apical','basal','perforant','soma']
stim_loc = ['apical','basal','perforant']
stim_dcs = ['control','cathodal','anodal']
drug = ['none']

# recording time before and after LTP induction
baset_pre = 20 
baset_post = 60

# organization of comments
com_columns= ['channel','block','sample','unknown','comment number']

#==============================================================================
# # preallocate
#==============================================================================
ind_block = np.empty([4,1])

#==============================================================================
# # loop over new files
#==============================================================================
for slice in dir_new:
    matfile = sio.loadmat(fpath_r+slice)
    
    # date
    date = int(slice[0:8])
    
    # sample rate
    fs = stats.mode(matfile['samplerate'],None).mode

    # length of each recording block
    l = matfile['dataend']-matfile['datastart']
    
    #==============================================================================
    # extract experiment info from comments
    #==============================================================================
    # loop over comments
    for idx,comment in enumerate(matfile['comtext']):
        com_chan = matfile['com'][idx,com_columns.index('channel')]
        com_block = matfile['com'][idx,com_columns.index('block')]        

        # recording location
        if 'Soma' in comment:
            chan_soma = com_chan
        if 'Apical' in comment:
            chan_apical = com_chan
        if 'Basal' in comment:
            chan_basal = com_chan
        if 'Perforant' in comment:
            chan_perforant = com_chan
            
        # induction blocks
        if 'induction_' in comment:
            ind_block[0] = com_block
        if 'induction2_' in comment:
            ind_block[1] = com_block
        if 'induction3_' in comment:
            ind_block[2] = com_block
        if 'induction4_' in comment:
            ind_block[3] = com_block
            

        # age
        if 'age' in comment:
            ageI = comment.index('age')
            ageL = len('age_')
            age = int(comment[ageI+ageL:ageI+ageL+2])
        
        # hemisphere
        if 'hemi' in comment:
            hemiI = comment.index('hemi')
            hemiL = len('hemi_')
            hemi = comment[hemiI+hemiL]
        
        # height
        if 'height' in comment:
            heightI = comment.index('height')
            heightL = len('height_')
            height = int(comment[heightI+heightL:heightI+heightL+2])
        
        # DCS current
        if 'current' in comment:
            currentI = comment.index('current')
            currentL = len('current_')
            current = int(comment[currentI+currentL:currentI+currentL+3])
        
    # bipolar pulse
    bipolar_date = 20161013 # changed the timing of bipolar pulse from 0 to 0.5 s
    if date >= bipolar_date:
        pulse_t = 0.5*fs
    else:
        pulse_t = 0
    
    base_index = ind_block[0]-baset_pre
            
            
            
            
        
    
    
    