clc
close all

% Directory where .mat files are located
mat_dir = '/Users/star/Documents/PNNL_Data/2024test'; 
 
% List both *nc.mat and *data.mat files
nc_mat_files = dir(fullfile(mat_dir, '*_nc.mat'));
data_mat_files = dir(fullfile(mat_dir, '*_data.mat'));

% Combine both file lists
mat_files = [nc_mat_files; data_mat_files];

% Initialize data storage
all_data = [];

displayDuration = 2; % Display Duration in seconds, adjust as needed

options.saveCSV = 0;  % Set to true to save CSV files
options.savePlots = 1; % Set to true to save plots
options.displayPlots = 1; % Set to true to display plots

%% -----Section 1-------------------Loop through each .mat file in the directory---------------------
start=1;
end_=288;

for k = 1:length(mat_files)
    mat_filename = fullfile(mat_files(k).folder, mat_files(k).name);
    MSave = load(mat_filename).MSave;
    
    % Extract the date from the filename
    [~, name, ~] = fileparts(mat_filename);
    
    % Find the date part in the filename using regular expressions
    date_match = regexp(name, '\d{8}', 'match');
    
    if ~isempty(date_match)
        date_str = date_match{1};
    else
        % Handle the case where a date is not found in the filename
        date_str = 'UnknownDate';
    end

    % Ensure required fields are present
    if ~isfield(MSave, 'nanosmps') || ~isfield(MSave.nanosmps, 'time')
        warning(['Data inconsistency in file: ', mat_filename]);
        continue;  % Skip this iteration
    end
   
% Extract times and interpolate variables
time_base = MSave.nanosmps.time / 3600;  % Convert to hours

% Check if time_base has 287 elements
if length(time_base) == 287
    % Prepend a 0 to make it 288 elements
    time_base = [0; time_base];
elseif length(time_base) == 289
    % If time_base has 289 elements, remove the last element to make it 288
    time_base = time_base(1:end-1);
elseif length(time_base) == 286
    % If time_base has 286 elements, add two elements to make it 288
    % add both at the beginning or both at the end, or one at each end
    % This example will add zeros at both ends
    time_base = [0; time_base; max(time_base)+1];  % max(time_base)+1 assumes uniform spacing and continues the pattern
end

    interpolated_vars = {
        'nanosmps', 'geometric_mean'; 
        'smps', 'geometric_mean';
        'nanosmps', 'total_N_conc';
        'smps', 'total_N_conc';
        'nanosmps', 'total_SA_conc';
        'smps', 'total_SA_conc';
        'ceilpblht', 'bl_height_1';
        'acsm', 'total_organics';
        'acsm', 'sulfate';
        'acsm', 'nitrate';
        'so2', 'so2';
        'met', 'rh_ambient';
        'met', 'temperature_ambient';
        'met', 'wind_direction';
        'met', 'wind_speed';
        'met', 'rain_amount';
        'ecorsf', 'turbulent_kinetic_energy';
        'ecorsf', 'co2_flux';
        'sirs20s', 'inst_diffuse';
        'qcrad1long','BestEstimate_down_short_hemisp';
    };

%% -----Section 2-----------------------storage all_data------------------------

% Pre-process the data
dN_dlogDp_nanosmps_interp = MSave.nanosmps.dN_dlogDp;
dN_dlogDp_nanosmps_interp(dN_dlogDp_nanosmps_interp == -9999) = NaN;  % Replace -9999 values with NaN

% Check the current number of columns
current_col_count = size(dN_dlogDp_nanosmps_interp, 2);

% Adjust column count based on the number present
if current_col_count == 289
    % Remove the last column if there are 289 columns
    dN_dlogDp_nanosmps_interp = dN_dlogDp_nanosmps_interp(:, 1:end-1);
elseif current_col_count == 287
    % If there are 287 columns, add a column of NaNs at the end
    dN_dlogDp_nanosmps_interp = [dN_dlogDp_nanosmps_interp, NaN(size(dN_dlogDp_nanosmps_interp, 1), 1)];
elseif current_col_count == 286
    % If there are 286 columns, add two columns of NaNs at the end
    dN_dlogDp_nanosmps_interp = [dN_dlogDp_nanosmps_interp, NaN(size(dN_dlogDp_nanosmps_interp, 1), 2)];
end

% Verify the size after adjustments
disp(['Adjusted size: ', num2str(size(dN_dlogDp_nanosmps_interp))]);

% Continue processing: Select only the desired rows (21 to 117) from dN_dlogDp_nanosmps_interp
selected_dN_dlogDp = dN_dlogDp_nanosmps_interp(21:117, :);


    % Initialize 'current_data'
    num_cols = size(interpolated_vars, 1) + 4; % +4: one for time, one for date, one for flag, and one for the smoothed data
    current_data = zeros(length(time_base), num_cols); 
    
    current_data(:, 1) = str2double(date_str);  % Insert the date column as the first column
    current_data(:, 2) = time_base;
    
    % Interpolate data
    for j = 1:size(interpolated_vars, 1)
        col_offset = 3;  % initial columns: date, time, flag
        if j > 1  
            col_offset = col_offset + 1; % account for the additional smoothed data column after the first variable
        end
        
        if isfield(MSave, interpolated_vars{j, 1}) && isfield(MSave.(interpolated_vars{j, 1}), interpolated_vars{j, 2})
            time_var = MSave.(interpolated_vars{j, 1}).time / 3600;  % Convert to hours
            data_var = MSave.(interpolated_vars{j, 1}).(interpolated_vars{j, 2});
            current_data(:, j + col_offset) = interp1(time_var, data_var, time_base, 'linear', NaN);
        else
            current_data(:, j + col_offset) = NaN;
        end
    end
    
    smoothed_geometric_mean_nanosmps = smooth(time_base, current_data(:, 4), 10, 'moving');
    current_data(:, 5) = smoothed_geometric_mean_nanosmps;
    
    % Calculate arithmetic mean for diameters below 30
    diameter_vals = MSave.nanosmps.diameter_mobility;
    valid_indices = diameter_vals < 15; % Indices of diameters < 15
    
    percentage_below_15 = zeros(length(time_base), 1); % Initialize the result array
    arithmetic_mean_below_30 = zeros(length(time_base), 1);
    
    for i = 1:length(time_base)
        % Get all concentrations at time i for diameters < 15 nm
        valid_concentrations = dN_dlogDp_nanosmps_interp(valid_indices, i);
        valid_diameters = diameter_vals(valid_indices);
    
        % Apply the mask to filter out NaN values
        valid_mask = ~isnan(valid_concentrations);
        product_sum = sum(valid_diameters(valid_mask) .* valid_concentrations(valid_mask));
    
        % Sum of concentrations for diameters < 15 nm at time i
        total_concentration_below_15 = sum(valid_concentrations(valid_mask));
        
        % Sum of all concentrations at time i (sum across all diameters)
        total_concentration = sum(dN_dlogDp_nanosmps_interp(:, i), 'omitnan'); % Summing all rows for column i
        
        % Check for non-zero total concentration before calculating the percentage
        if total_concentration ~= 0
            percentage_below_15(i) = (total_concentration_below_15 / total_concentration) * 100;
            arithmetic_mean_below_30(i) = product_sum / total_concentration_below_15;
    
        else
            percentage_below_15(i) = NaN; % If the total concentration is zero, the result is undefined
            arithmetic_mean_below_30(i) = NaN;
        end
    end

    % Now, percentage_below_15 is a time series of the percentage of particles below 15 nm.
    
    current_data(:, 17) = arithmetic_mean_below_30;
    current_data(:, 20) = percentage_below_15;

    % Compute the flag values
    smoothed_arithmetic_mean_below_30 = smooth(time_base, current_data(:, 17), 10, 'moving');
    current_data(:, 18) = smoothed_arithmetic_mean_below_30;
    smoothed_percentage_below_15 = smooth(time_base, current_data(:, 20), 10, 'moving');
    current_data(:, 21) = smoothed_percentage_below_15;
    
    % Call the updated function to calculate all three flags
    [flag1, flag2, flag3] = calculateFlags(current_data, percentage_below_15);
    [flag1, flag2, flag3] = postProcessFlags(flag1, flag2, flag3);
    
    % Assign the flag values to the appropriate columns in current_data
    current_data(:, 3) = flag1;  % Assign flag1 values
    current_data(:, 19) = flag2; % Assign flag2 values
    current_data(:, 22) = flag3; % Assign flag3 values
    
    % Set the shift amount (4 hours in terms of 5-min intervals)
    shift_amount = 4 * 60 / 5; % 4 hours shift
    
    % Initialize flag4 and flag5 with zeros
    flag4 = zeros(length(time_base), 1);
    flag5 = zeros(length(time_base), 1);
    
    % Perform the shift for flag4 (4 hours later)
    % Shift the entire flag3 column and fill the beginning with zeros
    flag4((shift_amount + 1):end) = flag3(1:(end - shift_amount));
    
    % Perform the shift for flag5 (4 hours earlier)
    % Shift the entire flag3 column and fill the end with zeros
    flag5(1:(end - shift_amount)) = flag3((shift_amount + 1):end);
    
    % Assign flag4 and flag5 values to the next available columns in current_data
    current_data(:, 23) = flag4; % Assign flag4 values
    current_data(:, 24) = flag5; % Assign flag5 values


 % Append the selected_dN_dlogDp data as new columns to current_data without transposing
    current_data(:, end+1:end+size(selected_dN_dlogDp', 2)) = selected_dN_dlogDp';

    % Initialize all variables with NaN to ensure they exist
bl_height_1_interp_full = NaN(length(time_base), 1);
total_organics_interp_full = NaN(length(time_base), 1);
sulfate_interp_full = NaN(length(time_base), 1);
nitrate_interp_full = NaN(length(time_base), 1);
rh_ambient_interp_full = NaN(length(time_base), 1);
temperature_ambient_interp_full = NaN(length(time_base), 1);
wind_direction_interp_full = NaN(length(time_base), 1);
wind_speed_interp_full = NaN(length(time_base), 1);
total_SA_conc_interp_full = NaN(length(time_base), 1);
so2_interp_full = NaN(length(time_base), 1);
turbulent_kinetic_energy_interp_full = NaN(length(time_base), 1);
BestEstimate_down_short_hemisp_interp_full = NaN(length(time_base), 1);

% Check and interpolate only if field exists
if isfield(MSave, 'ceilpblht')
    time_bl_height = MSave.ceilpblht.time / 3600;
    bl_height_1_interp_full = interp1(time_bl_height, MSave.ceilpblht.bl_height_1, time_base, 'linear', 'extrap');
end

if isfield(MSave, 'acsm')
    time_total_organics = MSave.acsm.time / 3600;
    total_organics_interp_full = interp1(time_total_organics, MSave.acsm.total_organics, time_base, 'linear', 'extrap');
    time_sulfate = MSave.acsm.time / 3600;
    sulfate_interp_full = interp1(time_sulfate, MSave.acsm.sulfate, time_base, 'linear', 'extrap');
    time_nitrate = MSave.acsm.time / 3600;
    nitrate_interp_full = interp1(time_nitrate, MSave.acsm.nitrate, time_base, 'linear', 'extrap');
end

if isfield(MSave, 'met')
    time_rh_ambient = MSave.met.time / 3600;
    rh_ambient_interp_full = interp1(time_rh_ambient, MSave.met.rh_ambient, time_base, 'linear', 'extrap');
    time_temperature_ambient = MSave.met.time / 3600;
    temperature_ambient_interp_full = interp1(time_temperature_ambient, MSave.met.temperature_ambient, time_base, 'linear', 'extrap');
    time_wind_direction = MSave.met.time / 3600;
    wind_direction_interp_full = interp1(time_wind_direction, MSave.met.wind_direction, time_base, 'linear', 'extrap');
    time_wind_speed = MSave.met.time / 3600;
    wind_speed_interp_full = interp1(time_wind_speed, MSave.met.wind_speed, time_base, 'linear', 'extrap');
end

if isfield(MSave, 'smps')
    time_total_SA_conc = MSave.smps.time / 3600;
    total_SA_conc_interp_full = interp1(time_total_SA_conc, MSave.smps.total_SA_conc, time_base, 'linear', 'extrap');
end

if isfield(MSave, 'so2')
    time_so2 = MSave.so2.time / 3600;
    so2_interp_full = interp1(time_so2, MSave.so2.so2, time_base, 'linear', 'extrap');
end

if isfield(MSave, 'ecorsf')
    time_turbulent_kinetic_energy = MSave.ecorsf.time / 3600;
    turbulent_kinetic_energy_interp_full = interp1(time_turbulent_kinetic_energy, MSave.ecorsf.turbulent_kinetic_energy, time_base, 'linear', 'extrap');
end

if isfield(MSave, 'qcrad1long')
    time_BestEstimate_down_short_hemisp = MSave.qcrad1long.time / 3600;
    BestEstimate_down_short_hemisp_interp_full = interp1(time_BestEstimate_down_short_hemisp, MSave.qcrad1long.BestEstimate_down_short_hemisp, time_base, 'linear', 'extrap');
end

    current_data(:, 9) = bl_height_1_interp_full;
    current_data(:, 10) = total_organics_interp_full;
    current_data(:, 11) = sulfate_interp_full;
    current_data(:, 12) = nitrate_interp_full;
    current_data(:, 13) = rh_ambient_interp_full;
    current_data(:, 14) = temperature_ambient_interp_full;
    current_data(:, 15) = wind_direction_interp_full;
    current_data(:, 16) = wind_speed_interp_full; 
%     current_data(:, 25) = total_SA_conc_nanosmps_interp;
    current_data(:, 26) = total_SA_conc_interp_full;
    current_data(:, 27) = so2_interp_full;
    current_data(:, 28) = turbulent_kinetic_energy_interp_full;
%     current_data(:, 29) = co2_flux_interp; 
%     current_data(:, 30) = inst_diffuse_interp; 
    current_data(:, 31) = BestEstimate_down_short_hemisp_interp_full; 
    

    % Append to the global data storage
    all_data = [all_data; current_data];
    % Replace -9999 values with NaN
    all_data(all_data == -9999) = NaN;
        
%% -----Section 3-----------------------------Common properties-----------------------------
    % Initialize time variables for each field to NaN
    time_hours_nanosmps = NaN(size(time_base));
    time_hours_smps = NaN(size(time_base));
%     time_hours_ccn = NaN(size(time_base));
    time_hours_acsm = NaN(size(time_base)); % Initialized to NaN
    time_hours_met = NaN(size(time_base));
    time_hours_ceilpblht = NaN(size(time_base));
    
    % Check and assign values if fields exist
    if isfield(MSave, 'nanosmps') && isfield(MSave.nanosmps, 'time')
        time_hours_nanosmps = MSave.nanosmps.time / 60 / 60;
    end
    if isfield(MSave, 'smps') && isfield(MSave.smps, 'time')
        time_hours_smps = MSave.smps.time / 60 / 60;
    end
    if isfield(MSave, 'ccn') && isfield(MSave.ccn, 'time_bounds')
    % Select the first row and transpose it to get a 110x1 vector
        time_column = MSave.ccn.time_bounds(1, :)';
    % Convert the time values to hours if they are in seconds
        time_hours_ccn = time_column / 60 / 60;
    end

    if isfield(MSave, 'acsm') && isfield(MSave.acsm, 'time')
        time_hours_acsm = MSave.acsm.time / 60 / 60; % This line is now protected by the if condition
    end
    if isfield(MSave, 'met') && isfield(MSave.met, 'time')
        time_hours_met = MSave.met.time / 60 / 60;
    end
    if isfield(MSave, 'so2') && isfield(MSave.so2, 'time')
        time_hours_so2 = MSave.so2.time / 60 / 60;
    end
    if isfield(MSave, 'ceilpblht') && isfield(MSave.ceilpblht, 'time')
        time_hours_ceilpblht = MSave.ceilpblht.time / 60 / 60;
    end
    if isfield(MSave, 'ecorsf') && isfield(MSave.ecorsf, 'time')
        time_hours_ecorsf = MSave.ecorsf.time / 60 / 60;
    end 
    if isfield(MSave, 'sirs20s') && isfield(MSave.sirs20s, 'time')
        time_hours_sirs20s = MSave.sirs20s.time / 60 / 60;
    end
    if isfield(MSave, 'qcrad1long') && isfield(MSave.qcrad1long, 'time')
        time_hours_qcrad1long = MSave.qcrad1long.time / 60 / 60;
    end

cycleLength = 5; % Assuming each cycle is of length 5
ccnDataExists = false; % Flag to indicate if ccn data exists

% Check if MSave has the 'ccn' field and the necessary subfields
if isfield(MSave, 'ccn') && all(isfield(MSave.ccn, {'supersaturation_calculated', 'N_CCN'}))
    ccnDataExists = true; % Set flag to true as ccn data exists

    % Determine the length to trim based on the length of time_column
    trimLength = length(time_column);

    % Trim MSave.ccn.supersaturation_calculated to match the length of time_column
    trimmedSupersaturation = MSave.ccn.supersaturation_calculated(1:trimLength);

    % Calculate the number of cycles
    numCycles = ceil(length(trimmedSupersaturation) / cycleLength);

    % Initialize arrays
    maxValues = zeros(1, numCycles);
    maxIndices = zeros(1, numCycles);
    correspondingTimes = zeros(1, numCycles);

    for i = 1:numCycles
        % Extract the current cycle
        startIdx = (i - 1) * cycleLength + 1;
        endIdx = min(i * cycleLength, length(trimmedSupersaturation));
        currentCycle = trimmedSupersaturation(startIdx:endIdx);

        % Find max value and index within this cycle
        [maxValues(i), idxInCycle] = max(currentCycle);
        maxIndices(i) = startIdx + idxInCycle - 1;

        % Extract the corresponding time
        correspondingTimes(i) = time_hours_ccn(maxIndices(i));
    end

    % Check if the indices are within bounds of N_CCN
    validIndices = maxIndices <= length(MSave.ccn.N_CCN);
    correspondingValues = NaN(1, numCycles);
    correspondingValues(validIndices) = MSave.ccn.N_CCN(maxIndices(validIndices));
else
    % If 'ccn' field does not exist, ensure the flag is set to false
    ccnDataExists = false;
    % Initialize correspondingValues and correspondingTimes to empty if needed
    correspondingTimes = [];
    correspondingValues = [];
end
        
    %% -----Section 4-----------------------------plot section-----------------------------------
            % Initial setup for consistent axis limits
        xlim_vals = [min(time_hours_nanosmps) max(time_hours_nanosmps)]; % Example x-limits, adjust as needed
        
        % Interpolating and smoothing data to align with time_hours_nanosmps
        % For SMPS
        if isfield(MSave, 'smps') && all(isfield(MSave.smps, {'dN_dlogDp', 'geometric_mean', 'total_N_conc', 'total_SA_conc'}))
            dN_dlogDp_smps_interp = interp1(time_hours_smps, MSave.smps.dN_dlogDp', time_hours_nanosmps)';
            geometric_mean_smps_interp = interp1(time_hours_smps, MSave.smps.geometric_mean, time_hours_nanosmps);
            total_N_conc_smps_interp = interp1(time_hours_smps, MSave.smps.total_N_conc, time_hours_nanosmps);
            total_SA_conc_smps_interp = interp1(time_hours_smps, MSave.smps.total_SA_conc, time_hours_nanosmps);
        else
            % Initialize to NaN if fields do not exist
            dN_dlogDp_smps_interp = NaN(size(time_hours_nanosmps));
            geometric_mean_smps_interp = NaN(size(time_hours_nanosmps));
            total_N_conc_smps_interp = NaN(size(time_hours_nanosmps));
            total_SA_conc_smps_interp = NaN(size(time_hours_nanosmps));
        end
        
        % For nanosmps
        if isfield(MSave, 'nanosmps') && all(isfield(MSave.nanosmps, {'geometric_mean', 'total_N_conc', 'total_SA_conc'}))
            geometric_mean_nanosmps_interp = interp1(time_hours_nanosmps, MSave.nanosmps.geometric_mean, time_hours_nanosmps);
            total_N_conc_nanosmps_interp = interp1(time_hours_nanosmps, MSave.nanosmps.total_N_conc, time_hours_nanosmps);
            total_SA_conc_nanosmps_interp = interp1(time_hours_nanosmps, MSave.nanosmps.total_SA_conc, time_hours_nanosmps);
        else
            % Initialize to NaN if fields do not exist
            geometric_mean_nanosmps_interp = NaN(size(time_hours_nanosmps));
            total_N_conc_nanosmps_interp = NaN(size(time_hours_nanosmps));
            total_SA_conc_nanosmps_interp = NaN(size(time_hours_nanosmps));
        end
        
        % For ceilpblht
        if isfield(MSave, 'ceilpblht') && isfield(MSave.ceilpblht, 'bl_height_1')
            bl_height_1_interp = interp1(time_hours_ceilpblht, MSave.ceilpblht.bl_height_1, time_hours_nanosmps);
        else
            bl_height_1_interp = NaN(size(time_hours_nanosmps));
        end
    
        % For so2
        if isfield(MSave, 'so2') && isfield(MSave.so2, 'so2')
            so2_interp = interp1(time_hours_so2, MSave.so2.so2, time_hours_nanosmps);
        else
            so2_interp = NaN(size(time_hours_nanosmps));
        end
        so2_interp_trimmed = so2_interp(6:end);
        time_hours_nanosmps_trimmed = time_hours_nanosmps (6:end);
    
    %     % For ccn
    %     if isfield(MSave, 'ccn') && isfield(MSave.ccn, 'N_CCN_dN')
    %         N_CCN_dN_interp = interp1(time_hours_ccn, MSave.ccn.N_CCN_dN', time_hours_nanosmps)';
    %         
    %     else
    %         N_CCN_dN_interp = NaN(size(time_hours_nanosmps));
    %     end
     if isfield(MSave, 'ccn') && isfield(MSave.ccn, 'N_CCN_dN') && isfield(MSave.ccn, 'time_bounds')
        % Determine the number of elements to extract from N_CCN based on the size of time_bounds
        num_elements = size(MSave.ccn.time_bounds, 2); 
    
        % Check if N_CCN has enough elements
        if length(MSave.ccn.N_CCN) >= num_elements
            % Extract the matching subset from N_CCN
            N_CCN_subset = MSave.ccn.N_CCN(1:num_elements); 
    
            % Interpolate N_CCN_dN
            N_CCN_dN_interp = interp1(time_hours_ccn, MSave.ccn.N_CCN_dN', time_hours_nanosmps)';
    
            % Interpolate the subset of N_CCN
            N_CCN_interp = interp1(time_hours_ccn, N_CCN_subset, time_hours_nanosmps);
        else
            error('MSave.ccn.N_CCN does not have enough data points for the desired subset.');
            N_CCN_dN_interp = NaN(size(time_hours_nanosmps));
            N_CCN_interp = NaN(size(time_hours_nanosmps));
        end
    else
        N_CCN_dN_interp = NaN(size(time_hours_nanosmps));
        N_CCN_interp = NaN(size(time_hours_nanosmps));
    end
    
    
        % For acsm
        if isfield(MSave, 'acsm') && all(isfield(MSave.acsm, {'total_organics', 'sulfate', 'nitrate'}))
            total_organics_interp = interp1(time_hours_acsm, MSave.acsm.total_organics, time_hours_nanosmps);
            sulfate_interp = interp1(time_hours_acsm, MSave.acsm.sulfate, time_hours_nanosmps);
            nitrate_interp = interp1(time_hours_acsm, MSave.acsm.nitrate, time_hours_nanosmps);
        else
            total_organics_interp = NaN(size(time_hours_nanosmps));
            sulfate_interp = NaN(size(time_hours_nanosmps));
            nitrate_interp = NaN(size(time_hours_nanosmps));
        end
        
        % For met
        if isfield(MSave, 'met') && all(isfield(MSave.met, {'rh_ambient', 'temperature_ambient', 'wind_direction', 'wind_speed','rain_amount'}))
            rh_ambient_interp = interp1(time_hours_met, MSave.met.rh_ambient, time_hours_nanosmps);
            temperature_ambient_interp = interp1(time_hours_met, MSave.met.temperature_ambient, time_hours_nanosmps);
            wind_direction_interp = interp1(time_hours_met, MSave.met.wind_direction, time_hours_nanosmps);
            wind_speed_interp = interp1(time_hours_met, MSave.met.wind_speed, time_hours_nanosmps);
            rain_amount_interp = interp1(time_hours_met, MSave.met.rain_amount, time_hours_nanosmps);
        else
            rh_ambient_interp = NaN(size(time_hours_nanosmps));
            temperature_ambient_interp = NaN(size(time_hours_nanosmps));
            wind_direction_interp = NaN(size(time_hours_nanosmps));
            wind_speed_interp = NaN(size(time_hours_nanosmps));
            rain_amount_interp = NaN(size(time_hours_nanosmps));
        end

        % For ecorsf
        if isfield(MSave, 'ecorsf') && all(isfield(MSave.ecorsf, {'turbulent_kinetic_energy', 'co2_flux'}))
            turbulent_kinetic_energy_interp = interp1(time_hours_ecorsf, MSave.ecorsf.turbulent_kinetic_energy, time_hours_nanosmps);
            co2_flux_interp = interp1(time_hours_ecorsf, MSave.ecorsf.co2_flux, time_hours_nanosmps);
        else
            turbulent_kinetic_energy_interp = NaN(size(time_hours_nanosmps));
            co2_flux_interp = NaN(size(time_hours_nanosmps));
        end
        
        % For sirs20s
        if isfield(MSave, 'sirs20s') && all(isfield(MSave.sirs20s, {'inst_diffuse'}))
            inst_diffuse_interp = interp1(time_hours_sirs20s, MSave.sirs20s.inst_diffuse, time_hours_nanosmps);  
        else
            inst_diffuse_interp = NaN(size(time_hours_nanosmps));
        end

         % For qcrad1long
        if isfield(MSave, 'qcrad1long') && all(isfield(MSave.qcrad1long, {'BestEstimate_down_short_hemisp'}))
            BestEstimate_down_short_hemisp_interp = interp1(time_hours_qcrad1long, MSave.qcrad1long.BestEstimate_down_short_hemisp, time_hours_nanosmps);  
        else
            BestEstimate_down_short_hemisp_interp = NaN(size(time_hours_nanosmps));
        end
 
 %%%%%%%%%%%%%-----------------------------------------------------------------------------------------------------------   
    if options.savePlots || options.displayPlots

        fig = figure('Position', [50, 50, 1800, 1000], 'Visible', options.displayPlots);
        set(fig, 'Renderer', 'opengl'); % or 'opengl'
        set(fig, 'GraphicsSmoothing', 'on'); % This is just an example

                global_font_size = 8; % Set font size only once
        set(gca, 'FontSize', global_font_size); % Set font size
    
                    % Adjust spacing between subplots
                    h = findall(gcf, 'Type', 'axes');  % This finds all axes objects in the figure
                    for i = 1:length(h)
                        pos = get(h(i), 'Position');
                        pos(2) = pos(2) + (i-1)*0.02; % Adjust the vertical spacing
                        set(h(i), 'Position', pos);
                    end
     
            % Adjust the end range based on the available data
            actual_end = min(end_ + 288*(k-1), size(all_data, 1));
   
         % Setup global plot aesthetics
set(groot, 'DefaultAxesFontSize', 12); % Set default font size for axes
set(groot, 'DefaultLineLineWidth', 1.5); % Set default line width   
    %-------------plot 1 Combined size distribution---------------
subplot(4,1,1);
% Plot H1 data
H1 = pcolor(time_hours_nanosmps, MSave.nanosmps.diameter_mobility, MSave.nanosmps.dN_dlogDp);
set(H1, 'EdgeColor', 'none');
hold on;

% Create a mask for overlay condition
overlay_mask = ~isnan(MSave.nanosmps.dN_dlogDp) & MSave.nanosmps.dN_dlogDp ~= -9999;

% Extend the mask to match the dimensions of H2 data
extended_mask = interp1(MSave.nanosmps.diameter_mobility, double(overlay_mask), MSave.smps.diameter_mobility, 'nearest', 'extrap');
extended_mask = logical(extended_mask);

% Adjust H2 data for overlay
dN_dlogDp_smps_interp(extended_mask) = NaN;

% Plot the adjusted H2 data
H2 = pcolor(time_hours_nanosmps, MSave.smps.diameter_mobility, dN_dlogDp_smps_interp);
set(H2, 'EdgeColor', 'none');

% Rest of the plotting details
caxis([0 4000]);
cbar = colorbar;
cbar.Label.String = {'d{\itN}/dlog_{10}{\itD}_p', '(cm^{-3})'};
ylabel('{\itD}_p (nm)', 'FontSize', 20);
yticks([5 20 100 400]);
set(gca, 'TickDir', 'out', 'Yscale', 'log', 'Xlim', [10 24], 'Ylim', [5 400], 'FontSize', 20);

box on;
set(gca, 'BoxStyle', 'full', 'XAxisLocation', 'bottom');

% Define x-axis tick positions and labels
xtick_positions = 10:2:24;  % Define tick positions from 10 to 24
xtick_labels = arrayfun(@(x) sprintf('%02d:00', x), xtick_positions, 'UniformOutput', false);  % Create labels formatted as HH:00
xticks(xtick_positions);
xticklabels(xtick_labels);

%-------------plot 2 Solar radiation and TKE---------------
subplot(4,1,2);
yyaxis left;
BestEstimate_down_short_hemisp_interp(BestEstimate_down_short_hemisp_interp == -9999) = NaN;
plot(time_hours_nanosmps, BestEstimate_down_short_hemisp_interp, 'DisplayName', 'SRI');
ylabel('SRI (W m^{-2})', 'FontSize', 20);

yyaxis right;
turbulent_kinetic_energy_interp(turbulent_kinetic_energy_interp == -9999) = NaN;
plot(time_hours_nanosmps, turbulent_kinetic_energy_interp, 'DisplayName', 'TKE', 'LineStyle', '--');
ylabel('TKE (m^2 s^{-2})', 'FontSize', 20);

xlim([10 24]);
legend('show', 'Location', 'eastoutside');
set(gca, 'FontSize', 20);

% Define x-axis tick positions and labels
xtick_positions = 10:2:24;  % Define tick positions from 10 to 24
xtick_labels = arrayfun(@(x) sprintf('%02d:00', x), xtick_positions, 'UniformOutput', false);  % Create labels formatted as HH:00
xticks(xtick_positions);
xticklabels(xtick_labels);

%-------------plot 3 RH and Temperature---------------
subplot(4,1,3);
rh_ambient_interp(rh_ambient_interp == -9999) = NaN;
temperature_ambient_interp(temperature_ambient_interp == -9999) = NaN;
yyaxis left;
plot(time_hours_nanosmps, rh_ambient_interp, 'DisplayName', 'RH');
ylabel('RH (%)', 'FontSize', 20);

yyaxis right;
plot(time_hours_nanosmps, temperature_ambient_interp, 'DisplayName', '{\it T}', 'LineStyle', '--');
ylabel('{\it T} (°C)', 'FontSize', 20);

xlim([10 24]);
legend('show', 'Location', 'eastoutside');
set(gca, 'FontSize', 20);

% Define x-axis tick positions and labels
xtick_positions = 10:2:24;  % Define tick positions from 10 to 24
xtick_labels = arrayfun(@(x) sprintf('%02d:00', x), xtick_positions, 'UniformOutput', false);  % Create labels formatted as HH:00
xticks(xtick_positions);
xticklabels(xtick_labels);

%-------------plot 4 Total Surface Area and Wind Direction---------------
subplot(4,1,4);
yyaxis left;
plot(time_hours_nanosmps, total_SA_conc_nanosmps_interp / 1000000, 'DisplayName', '{\it S}_{tot}');
ylabel('{\it S}_{tot} (\mum^2 cm^{-3})', 'FontSize', 20);

yyaxis right;
wind_direction_interp(wind_direction_interp == -9999) = NaN;
stem(time_hours_nanosmps, wind_direction_interp, 'filled', 'Marker', 'o', 'DisplayName', 'Wdir', 'LineStyle', 'none'); % Changed to stem plot
ylabel('Wdir', 'FontSize', 20);
ylim([0 360]);
yticks([0 90 180 270 360]);  % Set ticks at cardinal points
yticklabels({'N', 'E', 'S', 'W', 'N'});  % Label directions, note the sequence matches clockwise rotation

xlim([10 24]);
xlabel('Hour of Day (UTC)', 'FontSize', 20); % Set font size for x-axis label
legend('show', 'Location', 'eastoutside');
set(gca, 'FontSize', 20); % Adjust font size for all ticks

% Define x-axis tick positions and labels
xtick_positions = 10:2:24;  % Define tick positions from 10 to 24
xtick_labels = arrayfun(@(x) sprintf('%02d:00', x), xtick_positions, 'UniformOutput', false);  % Create labels formatted as HH:00
xticks(xtick_positions);
xticklabels(xtick_labels);

% Uniform alignment of subplot widths
align_width = get(subplot(4,1,2), 'Position');
for i = 1:4
    current_plot = subplot(4,1,i);
    pos = get(current_plot, 'Position');
    pos(3) = align_width(3) - 0.02; % Slightly reduce width for padding
    set(current_plot, 'Position', pos);
end

% Extract date from the .mat file name for use in the title and filename
date_str = regexp(mat_filename, '\d{8}', 'match', 'once');  % Extract YYYYMMDD pattern
% Adjust the title for the entire figure using the date_str
annotation('textbox', [0.5, 0.97, 0, 0], 'String', date_str, ...
           'FontSize', 25, 'FontWeight', 'bold', ...
           'EdgeColor', 'none', 'HorizontalAlignment', 'center');

    if options.savePlots
    % Save the plot using exportgraphics for higher quality output
            output_filename = fullfile(mat_dir, [date_str, '_Plot.png']);  % Save as PNG
            exportgraphics(gcf, output_filename, 'Resolution', 300); % Adjust the resolution as needed

    end
    drawnow; % Ensure the plot is rendered
        pause(displayDuration); % Wait for the specified duration
        close(fig); % Close the figure after the delay
    end

    % Clear large data variables at the end of each iteration to save memory
    clear MSave loaded_data;
end
%% ------section 5------------------------run flag----------------------

% Calculate start and end times based on flag values
start_time = NaN(size(all_data, 1), 1);
end_time = NaN(size(all_data, 1), 1);
for i = 2:size(all_data, 1)
    if all_data(i, 22) == 1 && all_data(i-1, 22) == 0
        start_time(i) = all_data(i, 2); % Mark start time
    elseif all_data(i, 22) == 0 && all_data(i-1, 22) == 1
        end_time(i-1) = all_data(i-1, 2); % Mark end time for the previous entry
    end
end

%%%%%%%----------------------------------------------------------------------------------------

if options.saveCSV
% Extracting the folder name from the folderPath
[~, folderName, ~] = fileparts(mat_dir);

% Generate the filename with the folder name
filename = fullfile(mat_dir , ['aligned_data_' folderName '.csv']);

% Rest for data processing
all_data = [all_data(:, 1:3), start_time, end_time, all_data(:, 4:end)];

% Array of particle sizes
particle_sizes = [2.0899999, 2.1700001, 2.2500000, 2.3299999, 2.4100001, 2.5000000, 2.5899999, 2.6900001, 2.7900000, 2.8900001, 3, 3.1099999, 3.2200000, 3.3399999, 3.4600000, 3.5899999, 3.7200000, 3.8499999, 4, 4.1399999, 4.2900000, 4.4499998, 4.6100001, 4.7800002, 4.9600000, 5.1399999, 5.3299999, 5.5200000, 5.7300000, 5.9400001, 6.1500001, 6.3800001, 6.6100001, 6.8499999, 7.0999999, 7.3699999, 7.6399999, 7.9099998, 8.1999998, 8.5100002, 8.8199997, 9.1400003, 9.4700003, 9.8199997, 10.200000, 10.600000, 10.900000, 11.300000, 11.800000, 12.200000, 12.600000, 13.100000, 13.600000, 14.100000, 14.600000, 15.100000, 15.700000, 16.299999, 16.799999, 17.500000, 18.100000, 18.799999, 19.500000, 20.200001, 20.900000, 21.700001, 22.500000, 23.299999, 24.100000, 25, 25.900000, 26.900000, 27.900000, 28.900000, 30, 31.100000, 32.200001, 33.400002, 34.599998, 35.900002, 37.200001, 38.500000, 40, 41.400002, 42.900002, 44.500000, 46.099998, 47.799999, 49.599998, 51.400002, 53.299999, 55.200001, 57.299999, 59.400002, 61.500000, 63.799999, 66.099998];

% Base headers before 'selected_dN_dlogDp' columns
header = {'date', 'time_hours', 'flag1', 'start_time', 'end_time', 'geometric_mean_nanosmps', 'smoothed_geometric_mean_nanosmps', 'geometric_mean_smps', 'total_N_conc_nanosmps', 'total_N_conc_smps', 'bl_height_1_interp', 'total_organics_interp', 'sulfate_interp', 'nitrate_interp', 'rh_ambient_interp', 'temperature_ambient_interp', 'wind_direction_interp', 'wind_speed_interp', 'arithmetic_mean_below_30', 'smoothed_arithmetic_mean_below_30', 'percentage below 15', 'smoothed percentage below 15', 'flag2', 'flag3', 'flag4', 'flag5'};

% Add dynamic headers for 'selected_dN_dlogDp' columns based on particle sizes
for i = 1:length(particle_sizes)
    % Format the particle size to have one decimal place and use it in the header
    header{end+1} = sprintf('Particle_Size_%.1fnm', particle_sizes(i));
end

% Now, 'header' includes dynamic headers for 'selected_dN_dlogDp' columns
csv_data = [header; num2cell(all_data)];

% Save the file with the new filename
writetable(cell2table(csv_data), filename, 'WriteVariableNames', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Helper function to calculate the flag value
% Given a matrix, it will iterate over its rows, computing the flag

function [flag1, flag2, flag3] = calculateFlags(all_data, percentage_below_15)
    % Initialize flags and state tracking variables
    numRows = size(all_data, 1);
    flag1 = zeros(numRows, 1);
    flag2 = zeros(numRows, 1);
    flag3 = zeros(numRows, 1);
    inFlag1State = false;
    inFlag2State = false;
    inFlag3State = false;

    for row = 2:numRows
        % Common calculations
        deltaT = all_data(row, 2) - all_data(row - 1, 2);
        rateChange = (deltaT ~= 0) * (all_data(row, 5) - all_data(row - 1, 5)) / deltaT;
        newRateChange = (deltaT ~= 0) * (all_data(row, 21) - all_data(row - 1, 21)) / deltaT;
        percentageCondition = percentage_below_15(row) > 10;

        % Calculations for flag1 and flag2
        for col = [5, 18] % Columns for flag1 and flag2
            deltaD = all_data(row, col) - all_data(row - 1, col);
            rateChange = (deltaT ~= 0) * deltaD / deltaT;

            % Start of flag1 or flag2
            if ~inFlag1State && all_data(row, col) < 9 && rateChange > 0
                if col == 5
                    flag1(row) = 1;
                    inFlag1State = true;
                elseif col == 18
                    flag2(row) = 1;
                    inFlag2State = true;
                end
            % Potential end of flag1 or flag2
            elseif (col == 5 && inFlag1State || col == 18 && inFlag2State) && rateChange <= 0
                % Check next 10 points
                futurePoints = min(row+10, numRows);
                futureRates = (all_data(row+1:futurePoints, col) - all_data(row:futurePoints-1, col)) ...
                            ./ (all_data(row+1:futurePoints, 2) - all_data(row:futurePoints-1, 2));
                
                % Determine whether to continue or end the flag
                if sum(futureRates > 0) >= 5
                    if col == 5
                        flag1(row) = 1;
                    elseif col == 18
                        flag2(row) = 1;
                    end
                else
                    if col == 5
                        inFlag1State = false;
                    elseif col == 18
                        inFlag2State = false;
                    end
                end
            % Continue flag1 or flag2
            elseif (col == 5 && inFlag1State || col == 18 && inFlag2State)
                if col == 5
                    flag1(row) = 1;
                elseif col == 18
                    flag2(row) = 1;
                end
            end
        end
      
       % Calculation for flag3
        if ~inFlag3State && all_data(row, 18) < 9 && newRateChange > 0
            flag3(row) = 1;
            inFlag3State = true; % Mark that flag3 has started
        elseif inFlag3State && rateChange <= 0
            % Check the next 10 points to determine if the flag should end
            futurePoints = min(row + 10, numRows);
            futureRates = (all_data(row+1:futurePoints, 18) - all_data(row:futurePoints-1, 18)) ...
                          ./ (all_data(row+1:futurePoints, 2) - all_data(row:futurePoints-1, 2));
            
            % If less than half of the next 10 rates are positive, end flag3
            if sum(futureRates > 0) < 5
                inFlag3State = false;
            else
                flag3(row) = 1; % Continue flag3
            end
        elseif inFlag3State
            flag3(row) = 1; % Continue flag3
        end
    end
        % Post-processing to fill gaps smaller than 12 time points
    gapStart = -1; % Initialize gap start index
    for row = 2:numRows
        if flag3(row) == 0 && flag3(row - 1) == 1
            gapStart = row; % Start of a gap
        elseif flag3(row) == 1 && gapStart ~= -1
            if row - gapStart < 12
                flag3(gapStart:row - 1) = 1; % Fill the gap
            end
            gapStart = -1; % Reset gap start index
        end
    end
end

    
%% ------section 6------------------------function----------------------
% Post-process to check for isolated flag=1 segments;
% this post-processing step is designed to remove short-lived segments where the flag is set to 1, 
% ensuring that any flagged segment in the final output has a length of at least 24 (two hour)

function [flag1, flag2, flag3] = postProcessFlags(flag1, flag2, flag3)
    % Adjust each flag array
    flag1 = adjustFlag(flag1);
    flag2 = adjustFlag(flag2);
    flag3 = adjustFlag(flag3);
end

function flag = adjustFlag(flag)
    count = 0;
    numRows = size(flag, 1);
    for row = 2:numRows
        if flag(row) == 1
            count = count + 1;
        else
            if count > 0 && count < 24 % Check if the previous segment was an isolated flag=1 segment
                flag(row - count : row - 1) = 0; % Reset the segment to flag=0
            end
            count = 0; % Reset count for the next segment
        end
    end
    % Check for the last segment in the array
    if count > 0 && count < 24
        flag(end - count + 1 : end) = 0;
    end
end

function alignLegendRight(lgd, rightmost_position, legend_width)
    lgd_pos = get(lgd, 'Position');
    lgd_pos(1) = rightmost_position - legend_width; 
    set(lgd, 'Position', lgd_pos);
end