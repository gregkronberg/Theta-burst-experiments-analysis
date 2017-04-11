# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 11:57:16 2017

@author: Greg Kronberg

Preprocessing raw data files for LTP + DCS experiments

Raw data are exported from LabChart as .mat files to the fpath_r directory

This script checks for new slices (against already processed slices in fpath_p directory)

Raw baseline data are organized as arrays with dimensions (time x recording blocks) 
(base_dend and base_soma for dendritic and somatic recordings, respectively)

Raw data during the plasticity induction period are stored separately (ind_dend and ind_soma)

Baseline traces are plotted for each slice and the window to take the slope is selected. 
Average slopes during this window are computed and normalized to the average of baseline slopes (stored as slopes_norm)
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
\\Theta LTP\\Processed Matlab Data'
elif comp==1: # desktop
    fpath_r = 'D:\\Google Drive\\Work\\Research Projects\\Theta LTP\
\\Raw Matlab Data\\'
    fpath_p = 'D:\\Google Drive\\Work\\Research Projects\\Theta LTP\
\\Processed Matlab Data\\'

# list files in each directory
dir_r = os.listdir(fpath_r)
dir_p = os.listdir(fpath_p)

# list new unprocessed data files
dir_new = list(set(dir_r)-set(dir_p))

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

# organization of comments
com_columns= ['channel','block','sample','unknown','comment number']

# preallocate
ind_block = np.empty([4,1],dtype=int)
loc_chan = -1*np.ones([len(rec_loc),1],dtype=int) # -1 means that location was not recorded
info = []
#==============================================================================
# loop over new files
#==============================================================================
for slice in dir_new:
    # load matlab file
    matfile = sio.loadmat(fpath_r+slice)
    
    # date
    date = slice[0:8]
    
    # sampling rate
    fs = stats.mode(matfile['samplerate'],axis=None).mode
                  
    # length of each recording block (samples)
    l_block = matfile['dataend']-matfile['datastart']
    
    #==========================================================================
    # extract experiment info from comments
    #==========================================================================
    # loop over comments
    for idx,comment in enumerate(matfile['comtext']):
        com_chan = matfile['com'][idx,com_columns.index('channel')]-1
        com_block = matfile['com'][idx,com_columns.index('block')]        
        
        # recording channels for each location
        for idx1,loc in enumerate(rec_loc):
            if loc in comment:
                loc_chan[idx1] = com_chan # 0 = channel 1, 1 = channel 2
        
        # recording block for each plasticity induction                
        for idx2,num in enumerate(ind_num):
            if num in comment:
                ind_block[idx2] = com_block
        
        # info about slice conditions
        for idx3,info_com in enumerate(slice_info): # loop over slice attributes
            if info_com in comment:
                info_start = re.search(info_com,comment).end()
                if info_com == 'height_':
                    info_end = info_start+3
                elif info_com == 'current_':
                    info_end = info_start+4
                else:
                    info_end1 = re.search(info_com+'.*?_',comment)
                    info_end = info_end1.end()
                info.append(comment[info_start:info_end-1]) # list of attributes (same order as slice_info)
                
    
    #==============================================================================
    # organize raw data into arrays
    # base_dend and base_soma are arrays of (time x blocks)
    #==============================================================================  
    # baseline index
    base_idx_pre = np.arange(ind_block[0]-baset_pre-1,ind_block[0]-1)
    base_idx_post = np.arange(ind_block[0],ind_block[0]+baset_post)
    base_idx = np.concatenate((base_idx_pre,base_idx_post),0)
    
    # organize baseline traces into column vectors
    # preallocate
    base_soma = np.empty([int(stats.mode(l_block,None).mode),
                          matfile['datastart'].shape[1]])
    base_dend = np.empty([int(stats.mode(l_block,None).mode),
                          matfile['datastart'].shape[1]])
    
    # loop over blocks
    for a in np.arange(0,matfile['datastart'].shape[1],dtype=int): 
        
        # check that it is a baseline block (not induction)
        if l_block[1,a] == stats.mode(l_block,None).mode: 
            
            # loop over possible recording locations
            for idx,loc in enumerate(rec_loc):
                # check for somatic recording
                if loc == 'Soma':
                    if loc_chan[idx] != -1:
                        # indeces for each block 
                        soma_idx = np.arange(matfile['datastart']
                            [int(loc_chan[idx]),a],matfile['dataend']
                            [int(loc_chan[idx]),a],dtype=int)
                        # store each block as a column
                        base_soma[:,a] = matfile['data'][:,soma_idx].T.reshape(-1) # (time x blocks)
                # repeat for dendritic recording
                elif loc_chan[idx] != -1:
                    dend_idx = np.arange(matfile['datastart'][loc_chan[idx],a],
                                     matfile['dataend'][loc_chan[idx],a],dtype=int)
                    base_dend[:,a] = matfile['data'][:,dend_idx].T.reshape(-1) # (time x blocks)
    
    # fill in gaps in recording
    # check for mismatch between number of recorded blocks and expected number
    if base_dend.shape[1]<max(base_idx):
        # repeat last block to fill gap
        base_dend = np.concatenate((base_dend,base_dend[:,-1].reshape(-1,1)),axis=1)
        base_idx = np.concatenate(base_idx,base_idx[-1]+1)
    if base_soma.shape[1]<max(base_idx):
        base_soma = np.concatenate((base_soma,base_soma[:,-1].reshape(-1,1)),axis=1)
        
    #==============================================================================
    # take slopes of baseline traces     
    #==============================================================================
    # bipolar pulse paramters
    pulse_date = 20161013 # changed timing of bipolar pulse after this date
    if date>= pulse_date:
        pulse_t = 0.5*fs
    else:
        pulse_t = 0
        
    # Time window to plot traces and select slopes
    slope_t1 = int(pulse_t + .001*fs)
    slope_t2 = int(pulse_t + .01*fs)
    
    # plot traces
    slope_plot = plt.figure()
    plt.plot(base_dend - base_dend[slope_t1,:])
    pylab.xlim([slope_t1,slope_t2])
    
    # select time window to take slope
    slope_win = plt.ginput(2)
    slope_win1 = int(slope_win[0][0])
    slope_win2 = int(slope_win[1][0])
    
    # take raw slope
    slopes_raw = np.mean(np.diff(base_dend[np.arange(slope_win1,slope_win2),:],
                                           axis=0),axis=0)*fs
    
    # normalize slopes
    slopes_base_mean = np.mean(slopes_raw[base_idx[0:baset_pre]])
    slopes_norm = slopes_raw/slopes_base_mean
    
    
    
        
        
        
                
            
        
    
    
    