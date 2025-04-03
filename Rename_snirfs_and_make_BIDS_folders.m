% Script that loops through each task in for a subject to rename snirf and
% move data to folders organized in BIDS format
% Assumes Replace_Probe.m is run first and new snrifs are created
%
close all
clear
clc

% cd '/projectnb/nphfnirs/ns/lcarlton/Homer3/'
% setpaths

%% Load file paths 
sub = 'BUMD006';    % CHANGE  
date = '2025-03-06'; % CHANGE

filePaths.PROJECTDIR = '/projectnb/nphfnirs/s/datasets/U01_ADRD/'; 
filePaths.RAWDATADIR = fullfile(filePaths.PROJECTDIR, 'sourcedata', 'raw');  % this is your RAW data.  
filePaths.SUBDIR = fullfile(filePaths.RAWDATADIR, sub, date);

run_folders = dir(filePaths.SUBDIR);
num_runs = length(run_folders) -2;    % count number of run folders (subtract 2 bc of how dir works)


% loop through each run folder and save updated snirf file
for i = 1:num_runs
    
    filePaths.RUNDIR = '';  % reset for each run
    filePaths.SNIRF = '';

    %run = char(strcat('00',num2str(i))); 
    date_run = run_folders(i+2).name;  % folder and file name format
    split_date_run = split(run_folders(i+2).name, '_');
    run = split_date_run{end};

    filePaths.RUNDIR = fullfile(filePaths.SUBDIR, date_run);    % path to folder for curr run
    filePaths.SNIRF = fullfile(filePaths.RUNDIR, char(strcat(date_run, '_SD.snirf')));   % curr snirf

    % load snirf
    if exist(filePaths.SNIRF, 'file') == 2
        disp(['Snirf file exists for run ', num2str(i)]);
    else
        disp(['Snirf file does not exist for run ', num2str(i)]);
        continue
    end
    snirfObj = SnirfLoad(filePaths.SNIRF);

    %% Create new BIDS compliant folder and copy snirf 
    % create new subject folder
    filePaths.BIDSDIR_SUB = fullfile( filePaths.PROJECTDIR, char(strcat('sub-', sub(end-1:end))), 'nirs');

    if ~exist(filePaths.BIDSDIR_SUB, 'dir')   % if directory does not exist, create it
        mkdir(filePaths.BIDSDIR_SUB)
        disp(['BIDS folder for sub-', sub(end-1:end), ' created'])
    end

    % create derivatives folder  if it doesn't already exist
    filePaths.DERIV = fullfile( filePaths.PROJECTDIR, 'derivatives');
    if ~exist(filePaths.DERIV, 'dir')   % if directory does not exist, create it
        mkdir(filePaths.DERIV)
        disp(['Derivatives folder created'])
    end

    %% Rename snirf based on description.json file
    % Load json file 
    filePaths.JSON = fullfile(filePaths.RUNDIR, char(strcat(date_run, '_description.json')));
    json_text = fileread(filePaths.JSON);
    description = jsondecode(json_text);

    json_split = split(description.experiment, '_');
    task_name = json_split{end};
    new_snirf_name = char(strcat('sub-', sub(end-1:end), '_task-', task_name, '_run-01_nirs.snirf'));
    
    filePaths.NEWSNIRF = fullfile(filePaths.BIDSDIR_SUB, new_snirf_name);
    
    % save snirf with BIDS formatting
    snirfObj.Save(filePaths.NEWSNIRF)

    %% Create table mapping each snirf to the original task number and description
    Subject = sub;
    Task_Number = date_run; Description_json = description.experiment;
    Snirf_file_name_BIDS = new_snirf_name;
    T{i} = table(string(Subject), string(Task_Number), string(Description_json), string(Snirf_file_name_BIDS));

end
%%
t = vertcat(T{:});
t.Properties.VariableNames = ["Subject", "Task_Number", "Description_json", "Snirf_file_name_BIDS"];
display(t)
writetable(t,fullfile(filePaths.DERIV, char(strcat('sub-', sub(end-1:end), '_snirf_task_info.txt'))),'Delimiter','\t')

