% Code to replace probe file bc its messed up for whatever reason
% before you run this code, you need to setpaths in homer3
%
% Loops through each task and replaces probe
%
% NOTE -> need to move all raw data to sourcedata folder to be BIDS
% compliant
% NOTE: code assumes task in description.json is under 'experiment' and is
% separated by '_'

% CODE NOTES:
% -make this 2 separate scripts
% -change derivatives folder where csv is saved
% -move creation of table to earlier in the code? and add stop point for if
%   its good and u continue or if u need to stop running to change something?
% -in new script, grab stims and make events.tsv, changing file names to be
%   meaningful
% -python script for DQR

% updated script 3/5/25
% Shannon Kelley
%clc
clear
close all


% cd '/projectnb/nphfnirs/ns/lcarlton/Homer3/'
% setpaths

%% Load file paths 
sub = 'BUMD003';    % CHANGE  
date = '2025-02-20'; % CHANGE

filePaths.PROJECTDIR = '/projectnb/nphfnirs/ns/Shannon/Dementia_project'; % CHANGE to your path (R:\KiranLab4\U01 Neuroscience of the Everyday World\8. fNIRS supplement ADRD\5. Data)
filePaths.RAWDATADIR = fullfile(filePaths.PROJECTDIR, 'RAW');  % this is your RAW data.   % CHANGE RAW to '2. fNIRS_data/'
filePaths.SUBDIR = fullfile(filePaths.RAWDATADIR, sub, date);
filePaths.PROBEDIR = '/projectnb/nphfnirs/ns/Shannon/Dementia_project/Code/probe_corrected.SD'; % CHANGE to your path

SDdata = load(filePaths.PROBEDIR,'-mat');   % load probe SD file

run_folders = dir(filePaths.SUBDIR);
num_runs = length(run_folders) -2;    % count number of run folders (subtract 2 bc of how dir works)

%% Create new directories for BIDS formatting 
filePaths.BIDSDIR = fullfile(filePaths.PROJECTDIR, 'Data');  % CHANGE 'Data' to something like '5. fNIRS_data_BIDS' 

if ~exist(filePaths.BIDSDIR, 'dir')   % if directory does not exist, create it
   mkdir(filePaths.BIDSDIR)
end


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
    
    %cd(filePaths.RUNDIR) %  cd into run folder

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
    
    %% Create new BIDS compliant folder and copy snirf 
    % create new subject folder
    filePaths.BIDSDIR_SUB = fullfile( filePaths.BIDSDIR, char(strcat('sub-', sub(end-1:end))), 'nirs');

    if ~exist(filePaths.BIDSDIR_SUB, 'dir')   % if directory does not exist, create it
        mkdir(filePaths.BIDSDIR_SUB)
        disp(['BIDS folder for sub-', sub(end-1:end), ' created'])
    end

    % create derivatives folder
    filePaths.DERIV = fullfile( filePaths.BIDSDIR, 'derivatives');
    if ~exist(filePaths.DERIV, 'dir')   % if directory does not exist, create it
        mkdir(filePaths.DERIV)
        disp(['Derivatives folder created'])
    end

    if ~exist(fullfile(filePaths.DERIV, char(strcat('sub-', sub(end-1:end)))), 'dir') % if directory does not exist, create it
        mkdir(fullfile(filePaths.DERIV, char(strcat('sub-', sub(end-1:end)))))
        disp(['Derivatives folder for sub-', sub(end-1:end), ' created'])
    end

    %copyfile(filePaths.SAVEPATH, filePaths.BIDSDIR_SUB)

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
t = vertcat(T{:});
writetable(t,fullfile(filePaths.DERIV, char(strcat('sub-', sub(end-1:end))), char(strcat('sub-', sub(end-1:end), 'snirf_task_info.csv'))),'Delimiter',' ')
