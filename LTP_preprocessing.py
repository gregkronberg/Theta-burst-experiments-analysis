# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 11:57:16 2017

@author: Greg Kronberg

Preprocessing raw data files for LTP + DCS experiments
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

#==============================================================================
# Set up directories
#==============================================================================
# choose computer to run on (0 = laptop, 1 = lab desktop)
comp = 1
if comp == 0: # laptop
    fpath_r = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data'
    fpath_p = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data'
elif comp==1: # desktop
    fpath_r = 'D:\\Google Drive\\Work\\Research Projects\\Theta LTP\\Raw Matlab Data\\'
    fpath_p = 'D:\\Google Drive\\Work\\Research Projects\\Theta LTP\\Processed Matlab Data\\'

# list files in each directory
dir_r = os.listdir(fpath_r)
dir_p = os.listdir(fpath_p)

# list new unprocessed data files
dir_new = list(set(dir_r)-set(dir_p))

#==============================================================================
# experimental conditions
#==============================================================================
rec_loc = ['apical','basal','perforant','soma']
stim_loc = ['apical','basal','perforant']
stim_dcs = ['control','cathodal','anodal']
drug = ['none']

# baseline recording period (minutes)
baset_pre = 20
baset_post = 60

# organization of comments
com_columns= ['channel','block','sample','unknown','comment number']

# preallocate
ind_block = np.empty([4,1],dtype=int)

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
    chan_soma=-1
    chan_apical=-1
    chan_basal=-1
    chan_perforant=-1
    for idx,comment in enumerate(matfile['comtext']):
        com_chan = matfile['com'][idx,com_columns.index('channel')]-1
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
        if 'age_' in  comment:
            ageI = comment.index('age_')
            ageL = len('age_')
            age = int(comment[ageI+ageL:ageI+ageL+2])
        
        # height of slice along dorsal/ventral axis
        if 'height_' in  comment:
            heightI = comment.index('height_')
            heightL = len('height_')
            height = int(comment[heightI+heightL:heightI+heightL+2])  
        
        # DCS current for 20 V/m field
        if 'current_' in  comment:
            currentI = comment.index('current_')
            currentL = len('current_')
            current = int(comment[currentI+currentL:currentI+currentL+3])

        # brain hemisphere
        if 'hemi_' in  comment:
            hemiI = comment.index('hemi_')
            hemiL = len('hemi_')
            hemi = comment[hemiI+hemiL:hemiI+hemiL+1]   
    
    #==============================================================================
    # organize raw data into arrays
    #==============================================================================
    # bipolar pulse paramters
    pulse_date = 20161013 # changed timing of bipolar pulse after this date
    if date>= pulse_date:
        pulse_t = 0.5*fs
    else:
        pulse_t = 0
        
    # baseline index
    base_idx_pre = np.arange(ind_block[0]-baset_pre,ind_block[0])
    base_idx_post = np.arange(ind_block[0]+1,ind_block[0]+baset_post+1)
    base_idx = np.concatenate((base_idx_pre,base_idx_post),0)
    
    # organize baseline traces into column vectors
    # preallocate
    base_soma = np.empty([int(stats.mode(l_block,None).mode),matfile['datastart'].shape[1]])
    base_dend = np.empty([int(stats.mode(l_block,None).mode),matfile['datastart'].shape[1]])
    
    # loop over blocks
    for a in np.arange(0,matfile['datastart'].shape[1],dtype=int): 
        
        # check that it is a baseline block (not induction)
        if l_block[1,a] == stats.mode(l_block,None).mode: 
            
            # check for somatic recording 
            if chan_soma != -1: 
                # indeces for each block 
                soma_idx = np.arange(matfile['datastart'][int(chan_soma),a],
                                     matfile['dataend'][int(chan_soma),a],dtype=int)
                # store each block as a column
                base_soma[:,a] = matfile['data'][:,soma_idx].T.reshape(-1)
            
            # repeat for dendritic recordings
            if chan_apical !=-1:
                dend_idx = np.arange(matfile['datastart'][chan_apical,a],
                                     matfile['dataend'][chan_apical,a],dtype=int)
                base_dend[:,a] = matfile['data'][:,dend_idx].T.reshape(-1)
            elif chan_basal !=-1:
                dend_idx = np.arange(matfile['datastart'][chan_basal,a],
                                     matfile['dataend'][chan_basal,a],dtype=int)
                base_dend[:,a] = matfile['data'][:,dend_idx].T.reshape(-1)
            elif chan_perforant !=-1:
                dend_idx = np.arange(matfile['datastart'][int(chan_perforant),a],
                                     matfile['dataend'][int(chan_perforant),a],dtype=int)
                base_dend[:,a] = matfile['data'][:,dend_idx].T.reshape(-1)
    
    # fill in gaps in recording
    if base_dend.shape[1]<max(base_idx):
        base_dend = np.concatenate((base_dend,base_dend[:,-1].reshape(-1,1)),axis=1)
    if base_soma.shape[1]<max(base_idx):
        base_soma = np.concatenate((base_soma,base_soma[:,-1].reshape(-1,1)),axis=1)
        
        
        
                
            
        
    
    
    