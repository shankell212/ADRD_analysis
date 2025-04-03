# -*- coding: utf-8 -*-
"""
Dementia Data analysis

Created on Mon Jan 20 11:55:31 2025

@author: shank
"""
#%% Import packages

import cedalion
import cedalion.nirs
import cedalion.sigproc.quality as quality
import cedalion.sigproc.motion_correct as motion_correct
import cedalion.xrutils as xrutils
import cedalion.datasets as datasets
import numpy as np
import xarray as xr
import pint
import matplotlib.pyplot as plt
import cedalion.plots as plots
from cedalion import units
import scipy.signal
import os.path
import pandas as pd
from cedalion.vis import plot_probe as vpp
from cedalion.vis import time_series as vts
#from circle_probe_cedalion import plot_circle_probe

import sys
sys.path.append('C:\\Users\\shank\\Documents\\GitHub\\IWHD_esplanade\\Python implentation\\DAB\\')
#import processingFuncs_DAB as pfDAB

sys.path.append('C:\\Users\\shank\\Documents\\GitHub\\cedalion-dab-funcs2')
import DABfuncs_group_avg as pfDAB_grp_avg
import DABfuncs_plot_DQR as pfDAB_dqr
import DABfuncs_load_and_preprocess as pfDABprocess


#%% User inputs
testing = 0
chan = "S10D94"

# Prune channels
snr_thresh = 10 
sd_threshs = [0, 45]*units.mm
amp_threshs = [0, 1e7] #*units.volts    

# Bandpass filter
fmin = 0.01 * units.Hz # change to 0 if doing GLM w/ drift correction
fmax = 0.5 * units.Hz

# SplineSG
flag_do_SplineSG = 0
frame_size = 10 * units.s
p = 0.99

# GLM with SSR
trange_hrf = [2, 25] * units.s # time range for block averaging
trange_hrf_stat = [5, 20] # time range for t-stat

stim_lst_hrf = [2,3] # for calculating HRFs

# this is utilized in the block average function if you do GLM analysis
ssr_rho_thresh = 15 * units.mm
glm_basis_func_param = 1 * units.s
glm_drift_order = 3
flag_do_GLM = True # if True, will do GLM analysis, otherwise does block average
                    # The GLM doesn't presently handle any NaN's it seems. But 
                    # pruned channels have NaN's.. curious.
                    
stim_lst_dqr = [2, 3] # for DQR plots

#%% Load in Data
subj = "01"
task = 'GEN1'

task_list = ["Cloudy", "GEN1", "GEN2", "PERS1CB", "PERS2CB", "PrestoCB"] # sub 1
#task_list = ["CloudyCB", "GEN1CB", "GEN2CB", "PERS1", "PERS2", "Presto"] # sub 2

run_num = 1
base_dir = 'D:\\fNIRS\\DATA\\Dementia\\'

print(f"Analyzing data for subject {subj}, run {run_num}, task {task}")

run = f"sub-{subj}_task-{task}_run-0{run_num}"
path2snirf = f"{base_dir}sub-{subj}\\{run}_nirs.snirf"  
path2events = f"{base_dir}sub-{subj}\\{run}_events.tsv"

filenm = f'sub-{subj}_task-{task}_nirs'

recordings = cedalion.io.read_snirf(path2snirf) # load in snirf\
events_df = pd.read_csv(path2events, sep='\t')

recTmp = recordings[0]

events_df.columns = recTmp.stim.columns # rename columns bc they dont match format needed in cedalion

recTmp.stim = events_df # update rec.stim with correct stim marks
stim= recTmp.stim

#%% Preprocess Data
# Prune channels
rec, chs_pruned, sci, psp = pfDABprocess.pruneChannels(recTmp, snr_thresh, sd_threshs, amp_threshs)

# convert to OD
rec = pfDABprocess.ODandGVTD(rec)

# Get the slope of 'od' before motion correction and any bandpass filtering
foo = rec['od'].copy()
foo = foo.pint.dequantify()
slope_base = foo.polyfit(dim='time', deg=1).sel(degree=1)
slope_base = slope_base.rename({"polyfit_coefficients": "slope"})
slope_base = slope_base.assign_coords(channel = rec['od'].channel)
slope_base = slope_base.assign_coords(wavelength = rec['od'].wavelength)

# Motion correct - splineSG
if flag_do_SplineSG:
    rec["od_splineSG"] = motion_correct.motion_correct_splineSG(rec["od"], p=0.99, frame_size=frame_size)
# Motion correct - TDDR
rec["od_tddr"] = motion_correct.tddr( rec["od"] )


# Get slopes after TDDR before bandpass filtering
slope_tddr = rec['od_tddr'].polyfit(dim='time', deg=1).sel(degree=1)
slope_tddr = slope_tddr.rename({"polyfit_coefficients": "slope"})
slope_tddr = slope_tddr.assign_coords(channel = rec['od_tddr'].channel)
slope_tddr = slope_tddr.assign_coords(wavelength = rec['od_tddr'].wavelength)

# GVTD for TDDR before bandpass filtering
amp_tddr = rec['od_tddr'].copy()
amp_tddr.values = np.exp(-amp_tddr.values)
rec.aux_ts['gvtd_tddr'], _ = quality.gvtd(amp_tddr)

# Plot DQR
pfDAB_dqr.plotDQR( rec = rec, chs_pruned = chs_pruned, slope = [slope_base, slope_tddr], filenm = filenm, flagSave = True, filepath = base_dir, stim_lst_str = stim_lst_dqr )


# bandpass filter od_tddr
rec['od_tddr_bpfilt'] = cedalion.sigproc.frequency.freq_filter(rec['od_tddr'], fmin, fmax)
if flag_do_SplineSG:
    rec['od_splineSG_bpfilt'] = cedalion.sigproc.frequency.freq_filter(rec['od_splineSG'], fmin, fmax)
    
# convert to conc
dpf = xr.DataArray(
    [1, 1],
    dims="wavelength",
    coords={"wavelength": rec["amp"].wavelength},
)

# calculate concentrations
rec["conc_tddr"] = cedalion.nirs.od2conc(rec["od_tddr_bpfilt"], rec.geo3d, dpf, spectrum="prahl")
if flag_do_SplineSG:
    rec["conc_splineSG"] = cedalion.nirs.od2conc(rec["od_splineSG_bpfilt"], rec.geo3d, dpf, spectrum="prahl")


#%% Block Avg

conc_epochs_tmp = rec["conc_tddr"].cd.to_epochs(
                            stim,  # stimulus dataframe
                            set(stim[stim.trial_type.isin(stim_lst_hrf)].trial_type), # select events
                            before=trange_hrf[0],  # seconds before stimulus
                            after=trange_hrf[1],  # seconds after stimulus
                        )
if flag_do_SplineSG:
    conc_epochs_tmp_spline = rec["conc_splineSG"].cd.to_epochs(
                                stim,  # stimulus dataframe
                                set(stim[stim.trial_type.isin(stim_lst_hrf)].trial_type), # select events
                                before=trange_hrf[0],  # seconds before stimulus
                                after=trange_hrf[1],  # seconds after stimulus
                            )

# BL subtract & block avg
# tddr
baseline_conc = conc_epochs_tmp.sel(reltime=(conc_epochs_tmp.reltime < 0)).mean('reltime')
conc_epochs = conc_epochs_tmp - baseline_conc
blockaverage = conc_epochs.groupby('trial_type').mean('epoch')

if flag_do_SplineSG:
    blockaverage_mean_tmp = blockaverage.assign_coords(trial_type=[x + '-tddr' for x in blockaverage.trial_type.values])
    blockaverage = xr.concat([blockaverage, blockaverage_mean_tmp],dim='trial_type')
    
    # spline
    baseline_conc = conc_epochs_tmp.sel(reltime=(conc_epochs_tmp_spline.reltime < 0)).mean('reltime')
    conc_epochs_spline = conc_epochs_tmp_spline - baseline_conc
    blockaverage = conc_epochs_spline.groupby('trial_type').mean('epoch')
    
    blockaverage_mean_tmp = blockaverage.assign_coords(trial_type=[x + '-splineSG' for x in blockaverage.trial_type.values])
    blockaverage = xr.concat([blockaverage, blockaverage_mean_tmp],dim='trial_type')


#%% Visualize
#vts.run_vis(rec)
#vpp.run_vis(blockaverage, rec.geo2d, rec.geo3d)

