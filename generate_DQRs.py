#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get DQR for individual subject

Created on Fri Mar 21 10:05:09 2025

@author: smkelley
"""

import os
import cedalion
import cedalion.nirs
import cedalion.sigproc.quality as quality
import cedalion.xrutils as xrutils


import xarray as xr
import matplotlib.pyplot as p
import cedalion.plots as plots
from cedalion import units
import numpy as np
import pandas as pd
from math import ceil

import sys
sys.path.append('/projectnb/nphfnirs/s/datasets/U01_ADRD/code/ADRD_analysis/modules')
import module_plot_DQR as dqr
import module_load_and_preprocess as preproc



def run_DQR_all_tasks(subj, task_lst, snr_thresh, sd_thresh, amp_thresh):   # !!! make task_lst non variable in future?
    '''
    Function that loops through each task in task_lst for a single subject and generates and saves
    the DQR plot for each task. 
    
    Inputs: 
        subj (string) : subject number (ex '01')
        task_lst (list) : of strings indicating each task name for a subject.
        snr_thresh (int) : threshold for signal to noise ratio
        sd_thresh (list) : list of int class of source detector separation thresholds. 
            first element lower thresh, second element upper threshold. ex [0, 45]. 
        amp_thresh (list) : list of int class of amplitude thresholds. 
        fmin (int) : lower threshold  value for bandpass filter
        fmax (int) : upper threshold value for bandpass filter
    
    '''
    
    root_dir = "/projectnb/nphfnirs/s/datasets/U01_ADRD/"
    stim_lst = ['ControlAnswer','ExpAnswer']  
    runNum = '01'

    
    cfg_prune = {
        'snr_thresh' : snr_thresh, # CHANGE
        'sd_threshs' : sd_thresh *units.mm, # CHANGE  # defines the lower and upper bounds for the source-detector separation that we would like to keep
        'amp_threshs' : amp_thresh, # CHANGE   # define whether a channel's amplitude is within a certain range
        'perc_time_clean_thresh' : 0.6, 
        'sci_threshold' : 0.6,
        'psp_threshold' : 0.1,
        'window_length' : 5 * units.s, 
        'flag_use_sci' : True,   
        'flag_use_psp' : False   
    }
    
    cfg_bandpass = { 
        'flag_bandpass_filter' : False,
        'fmin' : 0.01 * units.Hz, 
        'fmax' : 0.05 * units.Hz  
    }
            

    cfg_preprocess = {
        'flag_prune_channels' : True, 
        'median_filt' : 1, # CHANGE  # set to 1 if you don't want to do median filtering
        'cfg_prune' : cfg_prune, 
        'cfg_bandpass' : cfg_bandpass,
        'flag_do_GLM_filter' : False, 
    }
    
    # loop through each task
    for task in task_lst:
        
        # run DQR
        print(f'Generating DQR plot for task: {task}')
        run_DQR_single_task(root_dir, subj, task, runNum, stim_lst, cfg_prune, cfg_bandpass, cfg_preprocess)

            


def run_DQR_single_task(root_dir, subj, task, runNum, stim_lst, cfg_prune, cfg_bandpass, cfg_preprocess):
    '''
    Function to run DQR on a single task for a single subj. 
    
    '''
    
    subDir = os.path.join(root_dir, f'sub-{subj}', 'nirs')
    run_nm = f'sub-{subj}_task-{task}_run-{runNum}'
    
    snirf_path = os.path.join(subDir, run_nm + '_nirs.snirf' )
    
    # check if the snirf file exists
    if not os.path.exists( snirf_path ):
        print( f"Error: File {snirf_path} does not exist" )
    else:
        records = cedalion.io.read_snirf( snirf_path ) 
        rec = records[0]
    
    events_path = os.path.join(subDir, run_nm + '_events.tsv' )
    
    # check if the events.tsv file exists
    if not os.path.exists( events_path ):
        print( f"Error: File {events_path} does not exist" )
    else:
        stim_df = pd.read_csv(events_path, sep='\t' )
        rec.stim = stim_df
    
    #%% Preprocess 
    # create folder if it doesn't exist'
    der_dir = os.path.join(root_dir, 'derivatives', 'plots')
    if not os.path.exists(der_dir):
        os.makedirs(der_dir)
    der_dir = os.path.join(root_dir, 'derivatives', 'plots', 'DQR')
    if not os.path.exists(der_dir):
        os.makedirs(der_dir)
    
    # Preprocess data with median filt
    rec = preproc.preprocess( rec, cfg_preprocess['median_filt'] )
    
    # Prune channels
    rec, chs_pruned, sci, psp = preproc.pruneChannels( rec, cfg_preprocess['cfg_prune'] )
    pruned_chans = chs_pruned.where(chs_pruned != 0.4, drop=True).channel.values # get array of channels that were pruned
    
    
    # Calculate OD 
    # if flag pruned channels is True, then do rest of preprocessing on pruned amp, if not then do preprocessing on unpruned data
    if cfg_preprocess['flag_prune_channels']:
        rec["od"] = cedalion.nirs.int2od(rec['amp_pruned'])                
    else:
        rec["od"] = cedalion.nirs.int2od(rec['amp'])
        del rec.timeseries['amp_pruned']   # delete pruned amp from time series
    rec["od_corrected"] = rec["od"]    # need to reassign to new rec_str to work w/ code
    
    # Calculate GVTD on pruned data
    amp_masked = preproc.prune_mask_ts(rec['amp'], pruned_chans)  # use chs_pruned to get gvtd w/out pruned data (could also zscore in gvtd func)
    rec.aux_ts["gvtd"], _ = quality.gvtd(amp_masked) 
    
    
    lambda0 = amp_masked.wavelength[0].wavelength.values
    lambda1 = amp_masked.wavelength[1].wavelength.values
    snr0, _ = quality.snr(amp_masked.sel(wavelength=lambda0), cfg_preprocess['cfg_prune']['snr_thresh'])
    snr1, _ = quality.snr(amp_masked.sel(wavelength=lambda1), cfg_preprocess['cfg_prune']['snr_thresh'])
    
    
    dqr.plotDQR( rec, chs_pruned, cfg_preprocess, run_nm, root_dir, stim_lst )





    