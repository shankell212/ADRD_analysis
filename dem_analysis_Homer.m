%% Framework for analyzing esplanade experiment data
% 1. load the data 
% 2. preprocess the data
%   3.1. remove NaNs, Infs
%   3.2 prune channels
%   3.3 convert intensity to OD
%   3.4 splineSG for motion correction
%   3.5 convert OD to HbO/HbR
% 3. run GLM 
% 4. reorder output
% 5. save output -> to use in python 
% optionally - save the average HRF to be plotted directly in image space 

% change SD file for the diff systems -- cannot combine syst 1 and 2

% 0. add Homer to path and set preprocessing parameters
clear all 
clc

cd '/projectnb/nphfnirs/ns/lcarlton/Homer3/'
setpaths

%%
cfg.prune.SDrange = [0, 45]; % acceptable range of SD distances
cfg.prune.SNRthresh = 10;
cfg.prune.dRange = [0, 1e7]; % acceptable range of raw intensity values, amp = volts

cfg.spline.p = 0.99; % degree of spline function 
cfg.spline.FrameSize_sec = 10;
turnOn = 1;

cfg.bandpassfilt.fmin = 0;
cfg.bandpassfilt.fmax = 0.5;

cfg.glm.trange = [-2,25];       % time range for HRF 
cfg.glm.glmSolveMethod = 1;     % ordinary least squares method
cfg.glm.idxBasis = 1;           % type of basis func for HRD (1= sequence of Gaussian functions)
cfg.glm.paramsBasis = [1,1];    % params for basis function, specifies width of of Gauss and step btwn each
cfg.glm.rhoSD_ssThresh = 15;     % max distance for short sep regression (mm) -- WHY 20 AND NOT 10?
cfg.glm.flagNuisanceRMethod = 1; % uses most correlated SS channel for short sep regression
cfg.glm.driftOrder = 3;          % 3rd order polynomial for drift correction
cfg.glm.c_vector = 0;           

%% 1. load the data

subject_list = ["02"];

run_num = 1;
task = "CloudyCB";
%nTask = 2; % num tasks (or stims ?)

filePaths.baseDir_data = "/projectnb/nphfnirs/ns/Shannon/Dementia_project/Data/";

filePaths.probe = "/projectnb/nphfnirs/ns/Shannon/Dementia_project/Code/probe.SD";

for ss = 1:length(subject_list) 
    subj = subject_list(ss);
    fprintf('\n loading data for subject %s, task %s, run %s \n', subj, task, num2str(run_num))
    
    run = char(join(['sub-', subj, '_task-', task, '_run-0', num2str(run_num), '_nirs'],''));

    filePaths.subjDir = char(join([filePaths.baseDir_data, 'sub-', subj, '/'],''));
    
    filePaths.snirf = char(join([filePaths.subjDir, run, '.snirf'],''));

    filePaths.saveDir = char(join([filePaths.baseDir_data, 'derivatives/results/sub-', subj, '\'],''));
    if ~exist(filePaths.saveDir, 'dir')
       mkdir(filePaths.saveDir)
    end
    filePaths.outputSave_conc = char(join([filePaths.saveDir, run, '_concHRFs.mat' ],''));
    filePaths.outputSave_dodavg = char(join([filePaths.saveDir, run, '_dodavgimg.mat' ],''));
    filePaths.outputSave_concavg = char(join([filePaths.saveDir, run, '_concavgimg.mat' ],''));

    snirf = SnirfLoad(filePaths.snirf); % load the data
    %load(filePaths.probe, '-mat')
    
    %% 2. preprocessing 
    %%
    fprintf('run preprocessing \n')
    % remove NaN and Inf
    hmrR_PreprocessIntensity_NAN(snirf.data); % spline interp to replace any NaNs in data
    
    % prune channels (based on intensity and SD range)
    mlActAuto = hmrR_PruneChannels(snirf.data, snirf.probe, [], [], cfg.prune.dRange, cfg.prune.SNRthresh, cfg.prune.SDrange);
    
    % convert intensity to OD
    dod = hmrR_Intensity2OD(snirf.data); % dod = d data matrix is change in optical density
    dod.dataTimeSeries(isinf(dod.dataTimeSeries)) = 0; % get rid of inf vals
    
    % motion correction - splineSG
    dod = hmrR_MotionCorrectSplineSG(dod, mlActAuto, cfg.spline.p, cfg.spline.FrameSize_sec, turnOn);
    
    % filter the data - 0.01-0.5Hz
    dod = hmrR_BandpassFilt(dod, cfg.bandpassfilt.fmin, cfg.bandpassfilt.fmax);

    % convert OD to concentration 
    conc = hmrR_OD2Conc( dod, snirf.probe, [1 1] );
    
    
    %% 3. run GLM 
    fprintf('run GLM \n')
    [data_yavg, data_yavgstd, nTrials, data_ynew, data_yresid, data_ysum2, beta_blks, yR_blks, hmrstats] = hmrR_GLM(conc, snirf.stim, snirf.probe, mlActAuto, [], [], [], cfg.glm.trange, cfg.glm.glmSolveMethod, cfg.glm.idxBasis, ...
                                                                                                                    cfg.glm.paramsBasis, cfg.glm.rhoSD_ssThresh, cfg.glm.flagNuisanceRMethod, cfg.glm.driftOrder, cfg.glm.c_vector);
    filePaths.outputSave_data_yavg = char(join([filePaths.saveDir, run, '_data_yavg.mat' ],''));
    save(filePaths.outputSave_data_yavg, 'data_yavg')
    
    allSubj_datayavg.(join(['sub_', subj],'')) = data_yavg;

    %%
    % pass dataclasses to the plotprobe2
    data_yavg_plotprobe = data_yavg;
    cols_with_zeros = all(data_yavg_plotprobe.dataTimeSeries == 0);
    data_yavg_plotprobe.dataTimeSeries(:, cols_with_zeros) = NaN;
    
    % initialize data obj
    data_obj = snirf.copy;
    data_obj.data = data_yavg_plotprobe;
    
    % pass data_obj as an argument to PlotProbe2
    PlotProbe2(data_obj)

end

%% average data_yavg over subjects then call plotProbe2
%{
data_yavg_subjavg = data_yavg.copy();
data_yavg_subjavg.dataTimeSeries = zeros(size(data_yavg.dataTimeSeries));

for s = 1:length(subject_list)
    subj = subject_list(s);
    data_yavg = allSubj_datayavg.(join(['sub_', subj],'')); 


    data_yavg_subjavg.dataTimeSeries = data_yavg_subjavg.dataTimeSeries + data_yavg.dataTimeSeries;
    
    data_yavg_plotprobe = data_yavg;  


    %% individ subjects
    cols_with_zeros = all(data_yavg_plotprobe.dataTimeSeries == 0);
    data_yavg_plotprobe.dataTimeSeries(:, cols_with_zeros) = NaN;
    
    % initialize data obj
    data_obj = snirf.copy;
    data_obj.data = data_yavg_plotprobe;
    
    % pass data_obj as an argument to PlotProbe2
    PlotProbe2(data_obj)

end

data_yavg_subjavg.dataTimeSeries = data_yavg_subjavg.dataTimeSeries / length(subject_list);


%%
data_yavg_plotprobe = data_yavg_subjavg;        
cols_with_zeros = all(data_yavg_plotprobe.dataTimeSeries == 0);
data_yavg_plotprobe.dataTimeSeries(:, cols_with_zeros) = NaN;

% initialize data obj
data_obj = snirf.copy;
data_obj.data = data_yavg_plotprobe;

% pass data_obj as an argument to PlotProbe2
PlotProbe2(data_obj)
%}












