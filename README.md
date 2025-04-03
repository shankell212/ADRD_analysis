# ADRD_analysis
Scripts and functions to preprocess and analyze ADRD data. 

**files needed:**
- probe_corrected.SD   (the correct probe to be placed in snirf files)
- roi_coords_dem.csv   (file that maps channel values to different brain regions)

**scripts needed:**
- Replace_Probe.m
- Rename_snirfs_and_make_BIDS_folders.m
- rename_stims.m
- generate_DQRs.py   (python/ cedalion)
- analysis script (or can run analysis in Homer GUI)
*need Cedalion environment created


Preprocessing Scheme (updated for auto)
1  run Replace probe script   ( Replace_Probe.m )
2. run rename snirfs and create bids folders script
3. run rename stims script
4. generate DQR's (see below)
5. Run analysis in Homer GUI or script based

**Steps to generate DQRs:**
1) open a scc login node session
2) Run the following commands:
  2.2) module load miniconda
  2.3) conda activate cedalion
  2.4) cd /projectnb/nphfnirs/s/datasets/U01_ADRD/code/ADRD_analysis/
  
  2.5) python3 -c "from generate_DQRs import run_DQR_all_tasks; run_DQR_all_tasks( subj = '06', task_lst = ['Cloudy', 'GEN1', 'GEN2', 'PERS1', 'PERS2', 'Presto'], snr_thresh=10, sd_thresh= [0, 45], amp_thresh=[0, 1e7])"

NOTE: change the subject number, task_lst if it varies, snr_thresh, sd_thresh, and amp_thresh to values you want to use. 


Function description:
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
