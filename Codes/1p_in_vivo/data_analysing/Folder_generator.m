% Specify the path where folders will be created
FilePath = 'your_file_path_containing_the_image_series';

% Specify the source directory where the files are located
SourceDir = 'your_file_path_containing_the_image_series';  % Adjust this to the directory containing your .nd2 files

% Check if the specified path exists, if not, create it
if ~exist(FilePath, 'dir')
    mkdir(FilePath);
    fprintf('Directory "%s" created.\n', FilePath);
end

% Set the starting and ending numbers
start_num = 1;
end_num = 15;

% Set the number of digits for zero-padding (e.g., 3 digits for '001')
num_digits = 3;

% Generate and create the folders
for i = start_num:end_num
    % Create folder name with fixed number of digits and leading zeros
    folder_name = sprintf(['%0', num2str(num_digits), 'd'], i);
    
    % Full path of the new folder
    full_folder_path = fullfile(FilePath, folder_name);
    
    % Check if the folder already exists
    if ~exist(full_folder_path, 'dir')
        mkdir(full_folder_path);  % Create the folder
        fprintf('Folder "%s" created.\n', full_folder_path);
    else
        fprintf('Folder "%s" already exists.\n', full_folder_path);
    end
end

% Now, search for the files and move them into the corresponding folders

% Define the file pattern to search for .nd2 files
%filePattern = fullfile(SourceDir, '30s1.87camera*.nd2');
% Filename pattern: '30s1.87cameraXXX.nd2', 
% e.g.30s1.87camera001, should be imaged and present in order, 001,002,003 etc
filePattern = fullfile(SourceDir, 'your_file_prefix*.nd2');
fileList = dir(filePattern);

% Loop over each file
for k = 1:length(fileList)
    % Get the filename
    baseFileName = fileList(k).name;
    fullFileName = fullfile(SourceDir, baseFileName);
    
    % Extract XXX from filename using regular expressions
    % Filename pattern: '30s1.87cameraXXX.nd2'
    % XXX can be one or more digits
    tokens = regexp(baseFileName, '30s1\.87camera(\d+)\.nd2', 'tokens');
    % change according to your filename pattern
    if ~isempty(tokens)
        XXX = tokens{1}{1}; % Extracted number as string (e.g., '001', '1', etc.)
        
        % Pad XXX with leading zeros to match num_digits
        folderNum = sprintf(['%0', num2str(num_digits), 'd'], str2double(XXX));
        
        % Construct the target folder path
        targetFolder = fullfile(FilePath, folderNum);
        
        % Check if the target folder exists
        if exist(targetFolder, 'dir')
            % Move the file into the target folder
            movefile(fullFileName, targetFolder);
            fprintf('Moved file "%s" to folder "%s".\n', baseFileName, folderNum);
        else
            warning('Target folder "%s" does not exist for file "%s".', folderNum, baseFileName);
        end
    else
        warning('Filename "%s" does not match expected pattern.', baseFileName);
    end
end
