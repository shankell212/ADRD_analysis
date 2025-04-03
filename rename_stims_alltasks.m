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

sub = '06';    % CHANGE  

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


    else
        disp('Events file does not exist, creating events table to edit');
        
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
        events.trial_type = cell(length(trial_type), 1); % Preallocate as a cell array
        % Create the table
        events = table(onset, duration, amplitude, trial_type);
        
        % Display the table
        disp(events);
    
    end
    
    % if no stims, skip current task
    if isempty(events)
        disp(['No stims recorded for task ', files(j).name, ' Cannot create events.tsv.'])
        continue
    end
    
    %%
    
    orignames={'1','2'};  % in 1st 2 subs, 2=ExpAnswer 3=ControlAnswer
    
    % new names 
    newnames={'ControlAnswer','ExpAnswer'};   % 1 = control and 2 = exp for all subs 6+
    
    %specify the duration of each stim type
    lengths=[20, 20];  % CHECK IF THIS THE CURRECT DURATION 
    
    
    %% modify
    
    N=length(snirf.stim);
    
    %%
    n = length(events.trial_type);
    %events_table.trial_type = cell(height(events_table), 1); % Preallocate as a cell array
    
    % convert to cell
    if isnumeric(events.trial_type)
        events.trial_type = cellstr(num2str(events.trial_type)); % Convert numbers to cell array of strings
    end
    %
    for ki=1:n
        
        temp=events.trial_type;
        indice=find(contains(orignames,temp(ki)));
    
        %events_table.trial_type(ki) = [];
        events.trial_type(ki) = {newnames{indice}};
        
        if ismember('duration', events.Properties.VariableNames)
            events.duration(ki) = lengths(indice);
        else
            events.Duration(ki) = lengths(indice);
        end
    
    end
    %%
    if ismember('Trial_type', events.Properties.VariableNames)
        events.Trial_type = string(events.Trial_type);
    else
        events.trial_type = string(events.trial_type);
    end
    
    %% save data
    
    writetable(events, filePaths.events, 'FileType', 'text', 'Delimiter', '\t');
    disp('events.tsv file created successfully.');

end


