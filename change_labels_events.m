% Change stim names 
% Presto and Cloudy = 4 control, 4 tasks
% if CB that means first 4 tasks r control
% 3/6/35
% Shannon Kelley
clear
clc
close all

% cd '/projectnb/nphfnirs/ns/lcarlton/Homer3/' % CHANGE 
% setpaths

%% Set paths and load in data

sub = '04';    % CHANGE  

filePaths.baseDir = '/projectnb/nphfnirs/s/datasets/U01_ADRD/';  % CHANGE 

%run_name=['sub-', sub, '_task-', task, '_run-', run]; %file name 


filePaths.subDir = fullfile(filePaths.baseDir, ['sub-', sub], 'nirs');

cd(filePaths.subDir)
files = dir('*.snirf');  % Lists all .snirf files in the current folder
numFiles = length(files);  % Get the number of files

for j = 1:numFiles
    events = [];
    filePaths.snirf = fullfile(filePaths.subDir, files(j).name);  % full path 2 snirf
    filePaths.events = fullfile(filePaths.subDir, [files(j).name(1:end-11), '_events.tsv']); % full path 2 events
    
    % Load snirf file
    if exist(filePaths.snirf, 'file') == 2
        fprintf('\nSnirf file exists');
        snirf = SnirfLoad(filePaths.snirf);
    else
        fprintf('\nSnirf file does not exist');
    end
    
    %%
    disp(['Renaming stims for file ', files(j).name])

    % load events.tsv file
    if exist(filePaths.events, 'file') == 2
        fileInfo = dir(filePaths.events); 
        if fileInfo.bytes > 0
            disp('Events file exists.');
            events = readtable(filePaths.events, 'FileType', 'text', 'Delimiter', '\t');
        else
            disp('Events file exists but is empty, creating events table to edit');
            
            % Create events table to edit
            onset = [];  % Initialize empty arrays
            duration = [];
            amplitude = [];
            trial_type = [];
            
            % Loop through the stimuli
            for i = 1:length(snirf.stim)
                stim_data = snirf.stim(i).data; % Extract onset, duration, amplitude
                num_trials = size(stim_data, 1);
                
                % Append to the arrays
                onset = [onset; stim_data(:,1)];
                duration = [duration; stim_data(:,2)];
                amplitude = [amplitude; stim_data(:,3)];
                trial_type = [trial_type; repmat(str2double(snirf.stim(i).name), num_trials, 1)];
            end
            events = table(onset, duration, amplitude, trial_type);
            
            % Display the table
            disp(events);
        end

    end

    %%
    if ismember('Duration', events.Properties.VariableNames)
        events = renamevars(events, 'Duration', 'duration');
    end

    %%

    if ismember('Trial_type', events.Properties.VariableNames)
        events = renamevars(events, 'Trial_type', 'trial_type');
    end

    %%
    if ismember('Amplitude', events.Properties.VariableNames)
        events = renamevars(events, 'Amplitude', 'value');
    end

    %%
    if ismember('amplitude', events.Properties.VariableNames)
        events = renamevars(events, 'amplitude', 'value');
    end
    
    %%
    if ismember('Onset', events.Properties.VariableNames)
        events = renamevars(events, 'Onset', 'onset');
    end
    
    %% save data
    
    writetable(events, filePaths.events, 'FileType', 'text', 'Delimiter', '\t');
    disp('events.tsv file created successfully.');


end
