% Define directory where .nc and .cdf files are located
data_dir = '/Users/star/Documents/PNNL_Data/sirs20s'; 

% List both .nc and .cdf files
nc_files = dir(fullfile(data_dir, '*.nc'));
cdf_files = dir(fullfile(data_dir, '*.cdf'));

% Concatenate file lists
all_files = [nc_files; cdf_files];

% Data categories to extract
data_categories = {'acsm', 'met', 'nanosmps', 'smps', 'so2', 'ceilpblht', 'ccn', 'ecorsf', 'sirs20s','qcrad1long'};

% Loop through each file in the directory
for k = 1:length(all_files)
    file_path = fullfile(all_files(k).folder, all_files(k).name);
    
    % Extract the date from the file's name
    date_str = regexp(file_path, '\d{8}', 'match', 'once');
    if isempty(date_str)
        continue; % Skip this file if no valid date found
    end

    % Prepare the save file name
    save_filename = fullfile(data_dir, [date_str, '_data.mat']);

    % Initialize or Load MSave
    if isfile(save_filename)
        % Load the existing .mat file
        MSave = load(save_filename);
        MSave = MSave.MSave;
    else
        % Initialize MSave if the .mat file does not exist
        MSave = struct;
    end

    % Identify the data category from filename
    category = '';
    for c = 1:length(data_categories)
        if contains(file_path, data_categories{c})
            category = data_categories{c};
            break;
        end
    end
    if isempty(category)
        continue;  % Skip this file if it doesn't belong to the desired categories
    end
    
    % Read data from the file
    info = ncinfo(file_path);
    varNames = {info.Variables.Name};

for j = 1:length(varNames)
    newData = ncread(file_path, varNames{j});

    % Check if field exists and has compatible dimensions
    if isfield(MSave, category) && isfield(MSave.(category), varNames{j})
        existingData = MSave.(category).(varNames{j});
        if size(existingData, 2) == size(newData, 2)
            % Append new data to existing data if column dimensions match
            MSave.(category).(varNames{j}) = [existingData; newData];
        else
            % Handle incompatible dimensions here (e.g., skip or log a warning)
            fprintf('Warning: Incompatible dimensions for %s in %s\n', varNames{j}, file_path);
        end
    else
        % Create new field
        MSave.(category).(varNames{j}) = newData;
    end
end
    
    % Save the updated data to the .mat file
    save(save_filename, 'MSave');
end

