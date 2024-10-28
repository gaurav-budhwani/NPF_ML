clc;
clear;

% Base directory where the files are located
baseDir = '/Users/star/Downloads/';

% File prefix and file names
filePrefix = 'tdump.';
%fileSuffixes = {'00', '01', '02', '03', '04', '05', '06'};

fileSuffixes = {'16', '17', '18', '19', '20', '21', '22', '23'};

% Format specification for reading the files
formatSpec = '%6s%6s%6s%6s%6s%6s%6s%6s%8s%9s%9s%9s%9s%9s%s%[^\n\r]';

% Initialize containers for latitude, longitude, and height
Lat = {};
Lon = {};
Height = {};

% Read each file and concatenate the data
for k = 1:length(fileSuffixes)
    suffix = fileSuffixes{k};
    filename = [baseDir filePrefix suffix '.txt'];
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'ReturnOnError', false);
    fclose(fileID);
    
    % Convert strings to numeric values for Lat, Lon, and Height
    LatData = str2double(dataArray{10});
    LonData = str2double(dataArray{11});
    HeightData = str2double(dataArray{12});

    % Store data in cell arrays
    Lat{end+1} = LatData;
    Lon{end+1} = LonData;
    Height{end+1} = HeightData;
end

% Visualization
figure;
hold on;
colormap(parula); % Use a jet colormap
cmin = min(cellfun(@min, Height));
cmax = max(cellfun(@max, Height));
caxis([cmin, cmax]); % Set the color axis limits
myColormap = colormap; % Retrieve the current colormap

% Plot each file, segment by segment with interpolated colors based on height
for i = 1:length(fileSuffixes)
    for j = 1:length(Lat{i})-1
        % Calculate the average height for the current and next point
        avgHeight = mean([Height{i}(j), Height{i}(j+1)]);
        % Map the average height to the appropriate colormap index
        colorIndex = round(((avgHeight - cmin) / (cmax - cmin)) * size(myColormap, 1));
        colorIndex = max(colorIndex, 1); % Ensure index is at least 1
        colorIndex = min(colorIndex, size(myColormap, 1)); % Ensure index does not exceed colormap size
        plot(Lon{i}(j:j+1), Lat{i}(j:j+1), 'Color', myColormap(colorIndex, :), 'LineWidth', 4);
    end
    % Place a label at the midpoint of each trajectory
    midIndex = round(length(Lat{i}) / 2);
    text(Lon{i}(midIndex), Lat{i}(midIndex), [fileSuffixes{i}], 'FontSize', 25);
end

colorbar; % Display the color bar to indicate height variation
% Color bar settings
cbar = colorbar; % Display the color bar to indicate height variation
cbar.Label.String = 'Height (m AGL)'; % Labeling the color bar with height units
cbar.Label.FontSize = 37; % Increase font size for the color bar label
set(gca,'FontName','Helvetica','FontSize',37);

% Add the shapefiles for the USA and Mexico
shapefileDirs = {'/Users/star/Downloads/Shape_file_US_Mexico/USA', '/Users/star/Downloads/Shape_file_US_Mexico/Mexico'};
shapefileNames = {'s_11au16.shp', 'Mexico_Polygon.shp'};

for i = 1:length(shapefileDirs)
    cd(shapefileDirs{i});
    kml = shaperead(shapefileNames{i});
    for e = 1:length(kml)
        geoshow(kml(e).Y, kml(e).X, 'Color', 'k', 'LineWidth', 2);
    end
end

% Reset directory to base directory
cd(baseDir);

%title('Backward trajectories by 2018-12-20 with time and height variation','FontSize', 40);
title('Backward trajectories by 2020-05-18 with time and height variation','FontSize', 40);
xlabel('Longitude','FontSize', 40);
ylabel('Latitude','FontSize', 40);
axis equal;
axis([-130 -70 20 50]); % Setting consistent axis limits

hold off;
