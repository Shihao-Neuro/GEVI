% Define the directory you want to work with
dataDir = 'Your_file_path_containing_images_each_in_a_folder\Spike_Info'; % Replace with the path to your directory

% Define the prefix you want to add
prefix = '1214VS97_5_';  % Replace this with your desired prefix to add ahead to file e.g. spikeinfo_001  ->1214VS97_5_spikeinfo_001

% Get a list of all files matching 'spikeinfo_*.mat' in the specified directory
filePattern = fullfile(dataDir, 'spikeinfo_*.mat');
fileList = dir(filePattern);

% Loop over each found file and rename it
for k = 1:length(fileList)
    % Get the old file name and path
    oldName = fileList(k).name;
    oldFullPath = fullfile(dataDir, oldName);
    
    % Construct the new name and path
    newName = [prefix oldName];
    newFullPath = fullfile(dataDir, newName);
    
    % Rename the file
    movefile(oldFullPath, newFullPath);
    
    % Display a message
    fprintf('Renamed %s to %s\n', oldFullPath, newFullPath);
end
