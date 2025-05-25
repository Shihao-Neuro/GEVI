% Clean and clear all before run
clear all; close all; clc;

% Define the ROIs of interest
ROIs = [2];  % Modify this array to select different ROIs

% Get a list of all .mat files in the directory that contain 'spikeinfo_' in the file name
files = dir('*spikeinfo_*.mat');

% Check if any files are found
if isempty(files)
    error('No files matching the pattern ''*spikeinfo_*.mat'' were found in the current directory.');
end

% Sort the files to ensure they are in order
[~, idx] = sort({files.name});
files = files(idx);

% Extract the prefix from the first file name (everything before 'spikeinfo_')
prefix_parts = split(files(1).name, 'spikeinfo_');
prefix = prefix_parts{1};

% Number of sessions and number of ROIs
numSessions = length(files);
numROIs = length(ROIs);

% Initialize cell arrays to hold the extracted values
% Cell arrays of size N x Z
NoiseAmp1Hz = cell(numSessions, numROIs);
SNRmedian1Hz = cell(numSessions, numROIs);
SNRmean1Hz = cell(numSessions, numROIs);
SNRpeak1Hz = cell(numSessions, numROIs);
SpikeN = cell(numSessions, numROIs);

% Initialize cell arrays to hold concatenated SNRRawTrace1Hz
% One cell per ROI
SNRRawTrace1Hz_concat = cell(1, numROIs);

% Loop through each file (session)
for idx = 1:numSessions
    % Load the .mat file
    data = load(files(idx).name);
    
    % Check if 'Neuron' exists in the loaded data
    if isfield(data, 'Neuron')
        Neuron = data.Neuron;
        
        % Loop through each ROI of interest
        for roi_idx = 1:numROIs
            ROI = ROIs(roi_idx);
            
            % Check if Neuron has at least ROI elements
            if length(Neuron) >= ROI && isfield(Neuron(ROI), 'SpikeInfo')
                SpikeInfo = Neuron(ROI).SpikeInfo;
                
                % Extract NoiseAmpRaw1Hz
                if isfield(SpikeInfo, 'NoiseAmpRaw1Hz') && ~isempty(SpikeInfo.NoiseAmpRaw1Hz)
                    NoiseAmp1Hz{idx, roi_idx} = SpikeInfo.NoiseAmpRaw1Hz;
                else
                    NoiseAmp1Hz{idx, roi_idx} = 'NA';
                end
                
                % Extract SNRmedianRaw1Hz
                if isfield(SpikeInfo, 'SNRmedianRaw1Hz') && ~isempty(SpikeInfo.SNRmedianRaw1Hz)
                    SNRmedian1Hz{idx, roi_idx} = SpikeInfo.SNRmedianRaw1Hz;
                else
                    SNRmedian1Hz{idx, roi_idx} = 'NA';
                end
                
                % Extract SNRmeanRaw1Hz
                if isfield(SpikeInfo, 'SNRmeanRaw1Hz') && ~isempty(SpikeInfo.SNRmeanRaw1Hz)
                    SNRmean1Hz{idx, roi_idx} = SpikeInfo.SNRmeanRaw1Hz;
                else
                    SNRmean1Hz{idx, roi_idx} = 'NA';
                end
                
                % Extract SNRpeakRaw1Hz
                if isfield(SpikeInfo, 'SNRpeakRaw1Hz') && ~isempty(SpikeInfo.SNRpeakRaw1Hz)
                    SNRpeak1Hz{idx, roi_idx} = SpikeInfo.SNRpeakRaw1Hz;
                else
                    SNRpeak1Hz{idx, roi_idx} = 'NA';
                end
                
                % Extract SpikeN
                if isfield(SpikeInfo, 'SpikeN') && ~isempty(SpikeInfo.SpikeN)
                    SpikeN{idx, roi_idx} = SpikeInfo.SpikeN;
                else
                    SpikeN{idx, roi_idx} = 'NA';
                end
                
                % Extract SNRRawTrace1Hz
                if isfield(SpikeInfo, 'SNRRawTrace1Hz') && ~isempty(SpikeInfo.SNRRawTrace1Hz)
                    % Concatenate SNRRawTrace1Hz across sessions
                    if isempty(SNRRawTrace1Hz_concat{roi_idx})
                        SNRRawTrace1Hz_concat{roi_idx} = SpikeInfo.SNRRawTrace1Hz(:); % Ensure column vector
                    else
                        SNRRawTrace1Hz_concat{roi_idx} = [SNRRawTrace1Hz_concat{roi_idx}; SpikeInfo.SNRRawTrace1Hz(:)];
                    end
                end
            else
                % Neuron(ROI) does not exist or does not have SpikeInfo
                % Assign 'NA' to indicate missing data
                NoiseAmp1Hz{idx, roi_idx} = 'NA';
                SNRmedian1Hz{idx, roi_idx} = 'NA';
                SNRmean1Hz{idx, roi_idx} = 'NA';
                SNRpeak1Hz{idx, roi_idx} = 'NA';
                SpikeN{idx, roi_idx} = 'NA';
                % Do not modify SNRRawTrace1Hz_concat
            end
        end
    else
        % 'Neuron' not present in data
        % Assign 'NA' for all ROIs
        for roi_idx = 1:numROIs
            NoiseAmp1Hz{idx, roi_idx} = 'NA';
            SNRmedian1Hz{idx, roi_idx} = 'NA';
            SNRmean1Hz{idx, roi_idx} = 'NA';
            SNRpeak1Hz{idx, roi_idx} = 'NA';
            SpikeN{idx, roi_idx} = 'NA';
            % Do not modify SNRRawTrace1Hz_concat
        end
    end
end

% After extraction, you can convert concatenated traces to matrices if lengths match
% For this example, we keep them as cell arrays

% Save the variables to a .mat file
save([prefix 'ROIs_SNR_Features.mat'], 'NoiseAmp1Hz', 'SNRmedian1Hz', 'SNRmean1Hz', 'SNRpeak1Hz', 'SpikeN', 'SNRRawTrace1Hz_concat', 'ROIs');

% Display a message indicating completion
disp('Data extraction complete. Variables have been saved.');
