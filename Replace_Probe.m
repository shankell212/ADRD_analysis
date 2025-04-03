% Code to replace probe file bc its messed up for whatever reason
% before you run this code, you need to setpaths in homer3
%
% Loops through each task and replaces probe
%
% NOTE: code assumes task in description.json is under 'experiment' and is
% separated by '_'

% updated script 3/28/25
% Shannon Kelley
%clc
clear
close all


% cd '/projectnb/nphfnirs/ns/lcarlton/Homer3/'
% setpaths

%% Load file paths 
sub = 'BUMD006';    % CHANGE  
date = '2025-03-06'; % CHANGE

filePaths.PROJECTDIR = '/projectnb/nphfnirs/s/datasets/U01_ADRD/'; 
filePaths.RAWDATADIR = fullfile(filePaths.PROJECTDIR, 'sourcedata','raw');  % this is your RAW data.   %
filePaths.SUBDIR = fullfile(filePaths.RAWDATADIR, sub, date);
filePaths.PROBEDIR = fullfile(filePaths.PROJECTDIR, 'code','probe_corrected.SD'); % CHANGE to your path

SDdata = load(filePaths.PROBEDIR,'-mat');   % load probe SD file

run_folders = dir(filePaths.SUBDIR);
num_runs = length(run_folders) -2;    % count number of run folders (subtract 2 bc of how dir works)


%%
T = {};
% loop through each run folder and save updated snirf file
for i = 1:num_runs
    
    filePaths.RUNDIR = '';  % reset for each run
    filePaths.SNIRF = '';

    %run = char(strcat('00',num2str(i))); 
    date_run = run_folders(i+2).name;  % folder and file name format
    split_date_run = split(run_folders(i+2).name, '_');
    run = split_date_run{end};

    filePaths.RUNDIR = fullfile(filePaths.SUBDIR, date_run);    % path to folder for curr run
    filePaths.SNIRF = fullfile(filePaths.RUNDIR, char(strcat(date_run, '.snirf')));   % curr snirf
    
    %% Create table mapping each snirf to the original task number and description
    filePaths.JSON = fullfile(filePaths.RUNDIR, char(strcat(date_run, '_description.json')));
    json_text = fileread(filePaths.JSON);
    description = jsondecode(json_text);
    
    Subject = sub;
    Task_Number = date_run; Description_json = description.experiment;
    T{i} = table(string(Subject), string(Task_Number), string(Description_json));


    % load snirf
    if exist(filePaths.SNIRF, 'file') == 2
        disp(['Snirf file exists for run ', num2str(i)]);
    else
        disp(['Snirf file does not exist for run ', num2str(i)]);
        continue
    end
    snirfObj = SnirfLoad(filePaths.SNIRF);

    % Update positions and landmarks
    snirfObj.probe.sourcePos2D = SDdata.SD.SrcPos;
    snirfObj.probe.detectorPos2D = SDdata.SD.DetPos;
    snirfObj.probe.sourcePos3D = SDdata.SD.SrcPos3D;
    snirfObj.probe.detectorPos3D = SDdata.SD.DetPos3D;
    snirfObj.probe.landmarkPos2D = SDdata.SD.Landmarks2D.pos;
    snirfObj.probe.landmarkLabels = SDdata.SD.Landmarks.labels';    % note: labels dont match ?
    snirfObj.probe.landmarkPos3D = SDdata.SD.Landmarks.pos;

    % rename the snirf file name below to the name you want to save
    saveFile = char(strcat(date_run, '_SD.snirf'));
    filePaths.SAVEPATH = fullfile(filePaths.RUNDIR, saveFile);
    snirfObj.Save(filePaths.SAVEPATH)
    

end

%% Display table
t = vertcat(T{:});
t.Properties.VariableNames = ["Subject", "Task Number", "Description json"];
display(t)

%writetable(t,fullfile(filePaths.DERIV, char(strcat('sub-', sub(end-1:end))), char(strcat('sub-', sub(end-1:end), 'snirf_task_info.csv'))),'Delimiter',' ')
