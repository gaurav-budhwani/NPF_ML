% Define directory where .nc and .cdf files are located
data_dir = '/Users/star/Documents/PNNL_Data/add_0917'; 

% List both .nc and .cdf files
nc_files = dir(fullfile(data_dir, '*.nc'));
cdf_files = dir(fullfile(data_dir, '*.cdf'));

% Combine .nc and .cdf files into one array
all_files = [nc_files; cdf_files];

% Existing mat files
mat_files = dir(fullfile(data_dir, '*.mat'));

% Create a map to track existing .mat files by date
mat_file_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:length(mat_files)
    date_key = regexp(mat_files(i).name, '\d{8}', 'match');
    if ~isempty(date_key)
        mat_file_map(date_key{1}) = fullfile(mat_files(i).folder, mat_files(i).name);
    end
end

% Data categories to extract
data_categories = {'acsm', 'met', 'nanosmps', 'smps', 'so2', 'ceilpblht', 'ccn', 'ecorsf', 'sirs20s', 'qcrad1long'};

% Loop through each file in the directory
for k = 1:length(all_files)
    file_path = fullfile(all_files(k).folder, all_files(k).name);
    
    % Extract the date from the file's name
    date_str = regexp(file_path, '\d{8}', 'match', 'once');
    if isempty(date_str) || ~isKey(mat_file_map, date_str)
        continue;  % Skip files if no corresponding .mat file exists
    end

    % Load the existing .mat file
    save_filename = mat_file_map(date_str);
    MSave = load(save_filename);
    MSave = MSave.MSave;

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
            % Create new field if it doesn't exist in the struct
            MSave.(category).(varNames{j}) = newData;
        end
    end

    % Save the updated data to the .mat file
    save(save_filename, 'MSave');
end
