% Clear the workspace and close all figures
clear all;
close all;

% Add necessary paths to access subfunctions and example data
addpath("./subfunctions/");

% Set batch mode flag
IsBatchMode = true;

%% Define base file paths and parameters

BaseFilePath = 'Your_file_path_containing_images_each_in_a_folder\'; % Base path for images

%Create Spike_Info folder
FilePath2 = fullfile(BaseFilePath, 'Spike_Info\');
if ~exist(FilePath2, 'dir')
    mkdir(FilePath2);
end
% Define parameters for spike detection and image processing
Param.SampleRate = 638;             % Imaging frame rate in Hz
Param.SpikeTemplateLength = 3;      % Length of the spike template in time points
Param.SpikeTemplateN = 100;         % Number of spikes used for template generation
Param.CutOffFreq = 10;              % Cutoff frequency for high-pass filtering (Hz)
Param.SpikePolarity = -1;           % Spike polarity (-1 for negative spikes)
Param.MinSpikeTemplateN = [1 1 1 1]; % Minimum number of spikes to calculate template
Param.SNRList = [3.5 4 4.5 5];      % List of SNR thresholds for spike detection
Param.HeadTailSize = 200;           % Size of the head and tail in data trimming
Param.CellEnvSize = 15;             % Size extending from soma to estimate local background

% Load the Regions of Interest (ROIs) from a file (shared across all images)
ROIFileNamePrefix = [BaseFilePath '30s1.71camera1-52ROI']; % Shared ROI file name
load([ROIFileNamePrefix '.mat']); % Load ROIs

% Initialize structures to hold combined data across images
numROIs = length(ROIs);
numImages = 0; % To keep track of the number of images processed
AllSpikeTraces = {}; % Cell array to hold traces per ROI per image
AllSpikeIndices = {}; % Cell array to hold spike indices per ROI per image
AllImageLabels = {}; % Cell array to hold image labels

%% Define the list of image numbers to process
image_numbers = 1:52; % Modify this to include all the images you wish to process

% Loop over each image number
for idx = 1:length(image_numbers)
    XX_num = image_numbers(idx);
    % Format XX as a zero-padded three-digit string, e.g., '001', '002', etc.
    XX = sprintf('%03d', XX_num);
    XX_padded = XX; % 'XXX', matching the folder names
    fprintf('Processing image %s\n', XX_padded);
    
    % Define file paths and file name prefixes
    FilePath = [BaseFilePath  XX_padded '\']; % The path where the files are located
    FileNamePrefix = ['30s1.71camera' XX_padded]; % The prefix of the file name

    % Check if the directory exists
    if ~exist(FilePath, 'dir')
        fprintf('Directory %s does not exist. Skipping image %s.\n', FilePath, XX_padded);
        continue;
    end

    % Check if the required files exist
    pca_file = [FilePath  'pca_' XX_padded '.mat'];
    image_file = [FilePath FileNamePrefix '.mat'];

    if ~exist(pca_file, 'file')
        fprintf('File %s does not exist. Skipping image %s.\n', pca_file, XX_padded);
        continue;
    end

    if ~exist(image_file, 'file')
        fprintf('File %s does not exist. Skipping image %s.\n', image_file, XX_padded);
        continue;
    end

    % Load PCA noise components from a file
    Tmp = load(pca_file, 'NoisePCA');
    NoisePCA = Tmp.NoisePCA;
    
    % Load the corrected image stack
    load(image_file);
    
    % Remove noise components from the image stack using PCA
    ImStack = noise_pca_crt(ImStackCrt, NoisePCA);

    % Append the image label for later use in plots
    numImages = numImages + 1;
    AllImageLabels{numImages} = XX_padded;
    
    % Define file paths and file name prefixes
    FilePath = [BaseFilePath  XX_padded '\']; % The path where the files are located
    FileNamePrefix = ['30s1.71camera' XX_padded]; % The prefix of the file name
    
    % Load PCA noise components from a file
    Tmp = load([FilePath  'pca_' XX_padded '.mat'],'NoisePCA');
    NoisePCA = Tmp.NoisePCA;
    
    %%
    % Load the corrected image stack
    load([FilePath FileNamePrefix '.mat']);
    
    % Remove noise components from the image stack using PCA
    ImStack = noise_pca_crt(ImStackCrt, NoisePCA);
    
    % Add the mean image back to the denoised image stack to obtain the corrected images
    CrtImg = ones(size(ImStackCrt)) .* mean(ImStackCrt, 3) + double(ImStack);
    
    %%
    % Compute a preview image by calculating the standard deviation across frames
    ImPreview = std(single(ImStack(:,:,1:min(1000,size(ImStack,3)))), 0, 3);
    ImPreview = ImPreview / max(ImPreview(:)); % Normalize the preview image
    ImPreview = repmat(ImPreview, 1, 1, 3);    % Convert to RGB format for visualization
    
    % Initialize a mask covering all ROIs
    ROITot = ImPreview(:,:,1) * 0;
    for ii = 1:length(ROIs)
         ROITot(ROIs{ii}) = 1;
    end
    
    % Display the combined ROI mask (optional in batch mode)
    if ~IsBatchMode
        figure(1); imagesc(ROITot); axis image;
    end
    
    % Initialize an empty mask for individual ROIs
    ROIMask = ImPreview(:,:,1) * 0;
    
    % Initialize a structure array to hold neuron data
    Neuron = [];
    
    % Loop over each ROI to extract and analyze data
    for ii = 1:length(ROIs)
         Neuron(ii).ROI = ROIs{ii};    % Store the ROI indices for the neuron
         Neuron(ii).SpikeInfo = [];    % Initialize the spike information
         if ~isempty(Neuron(ii).ROI)
            disp(['roi # : ' num2str(ii)]);
    
            % Create a mask for the current ROI
            ROIMask = ROIMask * 0; 
            ROIMask(Neuron(ii).ROI) = 1;
    
            % Display the ROI overlayed on the preview image (optional in batch mode)
            if ~IsBatchMode
                ImPreview = display_roi(ImPreview, ROIMask);
                figure(1); imagesc(ImPreview); axis image; title(['roi # ' num2str(ii)]);
            end
    
            % Extract the data for the current ROI and background indices
            [DataROI, SomaIdx, BkgIdx] = extract_roi(Param, ImStack, ROIMask, ROITot);
    
            % Extract spike information from the ROI data
            [SpikeInfo, Trace, Bkg] = spike_extract(Param, DataROI, SomaIdx, BkgIdx, IsBatchMode);
            Neuron(ii).SpikeInfo = SpikeInfo;
            Neuron(ii).Trace = Trace;
            Neuron(ii).Bkg = Bkg;
    
            % Generate and save the plot for 'SNR-Normalized Raw Trace (1
            % Hz HPF) with spikes'  
            figure;
            plot(SpikeInfo.SNRRawTrace1Hz, 'Color', [0, 230/255, 0]);
            title(['ROI #' num2str(ii) ' - Image ' XX ' - SNR-Normalized Raw Trace (1 Hz HPF)']);
            xlabel('Frame');
            ylabel('SNR-Normalized Raw Trace (1 Hz HPF)');
            hold on;
            if ~isempty(SpikeInfo) && isfield(SpikeInfo, 'SpikeIdx')
                % Define x positions for spikes
                x_spikes = SpikeInfo.SpikeIdx;
                % Set constant y-value for strokes (e.g., 10% above max signal)
                y_constant = max(SpikeInfo.SNRRawTrace1Hz) * 1.1;
                % Define delta for the length of short strokes
                delta = (max(SpikeInfo.SNRRawTrace1Hz) - min(SpikeInfo.SNRRawTrace1Hz)) * 0.04;
                y_start = y_constant - delta;
                y_end = y_constant + delta;
                % Create arrays for plotting
                x_positions = [x_spikes'; x_spikes'];
                y_positions = [y_start * ones(1, length(x_spikes)); y_end * ones(1, length(x_spikes))];
                % Plot short vertical strokes in darker blue color
                plot(x_positions, y_positions, 'Color', [0.3, 0.5, 1], 'LineWidth', 0.5);
            end
            hold off;
            % Save the figure
            saveas(gcf, [FilePath 'ROI_' num2str(ii) '_Trace_Spikes_Image_' XX '.png']);
            close(gcf); % Close the figure to prevent too many open figures
    
            % Collect traces and spike indices for later plotting
            AllSpikeTraces{ii, numImages} = SpikeInfo.SNRRawTrace1Hz;
            AllSpikeIndices{ii, numImages} = SpikeInfo.SpikeIdx;
    
        end
    end
    
    % Count the number of neurons with detected spikes
    NumberSpikeNeuron = count_spike_neuron(Neuron);
    
    % Save the neuron data and parameters to a file
    save([FilePath2  'spikeinfo_' XX_padded '.mat'], 'Neuron', 'Param');
    
end   % End of image_numbers loop

%% After processing all images, generate the plots for each ROI

for ii = 1:numROIs
    figure;
    numSubplots = numImages; % Number of images processed
    % Determine the subplot layout (e.g., 2x4 for 8 images)
    numRows = ceil(sqrt(numSubplots));
    numCols = ceil(numSubplots / numRows);
    
    for idx = 1:numSubplots
        subplot(numRows, numCols, idx);
        if size(AllSpikeTraces, 1) >= ii && size(AllSpikeTraces, 2) >= idx
            Trace = AllSpikeTraces{ii, idx};
            SpikeIdx = AllSpikeIndices{ii, idx};
        else
            Trace = [];
            SpikeIdx = [];
        end
        if isempty(Trace)
            % No data for this ROI in this image
            title(['Image ' AllImageLabels{idx} ' - No Data']);
            continue;
        end
        %plot(Trace, 'Color', [0, 230/255, 0]); % For VS74
        plot(Trace, 'Color', [255/255, 169/255, 64/255]); % For VS91
        hold on;
        if ~isempty(SpikeIdx)
            y_constant = max(Trace) * 1.1;
            delta = (max(Trace) - min(Trace)) * 0.04;
            y_start = y_constant - delta;
            y_end = y_constant + delta;
            x_positions = [SpikeIdx'; SpikeIdx'];
            y_positions = [y_start * ones(1, length(SpikeIdx)); y_end * ones(1, length(SpikeIdx))];
            plot(x_positions, y_positions, 'Color', [0.3, 0.5, 1], 'LineWidth', 0.5);
        end
        xlabel('Frame');
        ylabel('SNR-Normalized Raw Trace (1 Hz HPF)');
        title(['Image ' AllImageLabels{idx}]);
        ylim([-5, 10]);
        hold off;
    end
    % Replace suptitle with sgtitle or annotation
    % Option 1: Use sgtitle if available
    % sgtitle(['ROI #' num2str(ii)]);
    % Option 2: Use annotation if sgtitle is not available
    % Add a super title above all subplots
    annotation('textbox', [0 0.95 1 0.05], ...
        'String', ['ROI #' num2str(ii)], ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 12, 'FontWeight', 'bold');
    % Save the figure
    saveas(gcf, [BaseFilePath 'ROI_' num2str(ii) '_AllImages_Traces.png']);
    %close(gcf); % Close the figure to prevent too many open figures
end

%% Function to count neurons with detected spikes

function NumberSpikeNeuron = count_spike_neuron(Neuron)
    Count = 0;
    for ii = 1:size(Neuron, 2)
        if ~isempty(Neuron(ii).SpikeInfo)
            Count = Count + 1;
        end
    end
    NumberSpikeNeuron = Count;
end

%% Function to remove noise components using PCA correction

function ImgCrt = noise_pca_crt(ImgRaw, NoisePCA)
    NoiseCrtN = 1; % Number of PCA components to remove
    Trace = NoisePCA.Trace(1:NoiseCrtN, :); % Select the first few components

    ImgSize = size(ImgRaw);
    ImgRaw = reshape(single(ImgRaw), ImgSize(1)*ImgSize(2), ImgSize(3)); % Flatten spatial dimensions
    Coef = inv(Trace*Trace') * (Trace * ImgRaw'); % Calculate coefficients for noise components
    ImgCrt = ImgRaw - Coef' * Trace; % Subtract noise components from the images
    ImgCrt = uint16(ImgCrt - min(ImgCrt(:))); % Adjust pixel values to be positive integers
    ImgCrt = reshape(ImgCrt, ImgSize(1), ImgSize(2), ImgSize(3)); % Reshape back to original dimensions
end

%% Function to display the ROI overlayed on the preview image

function Img = display_roi(Img, ROIMask)
    ROIMask = repmat(ROIMask, 1, 1, 3);        % Convert ROI mask to RGB format
    ROIMask = ROIMask .* rand(1, 1, 3) / 6;    % Assign random color to ROI
    Img = Img + ROIMask;                       % Overlay ROI on the image
end

%% Function to extract ROI data, soma indices, and background indices

function [DataROI1, SomaIdx, BkgIdx] = extract_roi(Param, ImStack1, ROIMask, ROITot)
    SpikePolarity = Param.SpikePolarity;   % Polarity of spikes
    CellEnvSize = Param.CellEnvSize;       % Size extending from soma for background estimation

    % Find the bounding box of the ROI
    Tmp = find(max(ROIMask, [], 1));
    xmin = min(Tmp);
    xmax = max(Tmp);
    Tmp = find(max(ROIMask, [], 2));
    ymin = min(Tmp);
    ymax = max(Tmp);

    % Expand the bounding box by CellEnvSize, ensuring it stays within image bounds
    xmin = max(1, xmin - CellEnvSize);
    xmax = min(size(ROIMask, 2), xmax + CellEnvSize);
    ymin = max(1, ymin - CellEnvSize);
    ymax = min(size(ROIMask, 1), ymax + CellEnvSize);

    % Extract the image data within the extended bounding box
    DataROI1 = double(ImStack1(ymin:ymax, xmin:xmax, :));

    % Create masks for the soma and background within the bounding box
    SomaMask = ROIMask(ymin:ymax, xmin:xmax);
    BkgMask = 1 - ROITot(ymin:ymax, xmin:xmax);

    % Multiply data by spike polarity to correct for negative spikes
    DataROI1 = DataROI1 * SpikePolarity;

    % Find the indices of the soma and background pixels
    SomaIdx = find(SomaMask(:));
    BkgIdx = find(BkgMask(:));

    % Display the soma and background masks (optional in batch mode)
    if ~evalin('caller', 'IsBatchMode')
        figure(5000);
        subplot(1, 2, 1); imagesc(SomaMask); axis image; title('Soma Mask');
        subplot(1, 2, 2); imagesc(BkgMask); axis image; title('Background Mask');
    end
end

%% Function to extract spike information from ROI data

function [SpikeInfo, Trace, Bkg] = spike_extract(Param, DataROI1, SomaIdx, BkgIdx, IsBatchMode)

    % Create a mask for the soma
    SomaMask = DataROI1(:,:,1) * 0;
    SomaMask(SomaIdx) = 1;

    % Apply high-pass filtering at 25 Hz to the data (existing code)
    DataROIFilt = highpass_video_filt(Param, DataROI1);
    DataROIFilt = reshape(DataROIFilt, [], size(DataROIFilt, 3));

    % Compute the mean trace for the soma and background (existing code)
    RawTrace = mean(DataROIFilt(SomaIdx, :), 1);
    RawTrace = RawTrace(:);
    Bkg = mean(DataROIFilt(BkgIdx, :), 1);
    Bkg = Bkg(:);

    % Subtract background from the soma trace (existing code)
    Trace = RawTrace - Bkg;
    Trace = Trace - mean(Trace(:));

    % Estimate noise amplitude (existing code)
    Mask = (Trace < 0);
    NoiseAmp1 = sqrt(sum(Trace.^2 .* Mask) / sum(Mask));

    % Compute SNR-normalized raw trace (existing code)
    SNRRawTrace = Trace / NoiseAmp1;

    % Perform spike denoising and detection (existing code)
    SpikeInfo = spike_denoise(Param, Trace);

    % Compute SNR-normalized filtered trace (existing code)
    if ~isempty(SpikeInfo) && isfield(SpikeInfo, 'SNRFiltTrace')
        SNRFiltTrace = SpikeInfo.SNRFiltTrace;
    else
        % If no spikes detected, set SNRFiltTrace to SNRRawTrace
        SNRFiltTrace = SNRRawTrace;
    end

    % -------------------- Existing Code Above --------------------

    % Create a modified Param structure for 1 Hz high-pass filtering
    Param1Hz = Param;
    Param1Hz.CutOffFreq = 1; % Set cutoff frequency to 1 Hz

    % Apply high-pass filtering at 1 Hz to the original data
    DataROIFilt1Hz = highpass_video_filt(Param1Hz, DataROI1);
    DataROIFilt1Hz = reshape(DataROIFilt1Hz, [], size(DataROIFilt1Hz, 3));

    % Compute the mean trace for the soma and background (1 Hz filtered data)
    RawTrace1Hz = mean(DataROIFilt1Hz(SomaIdx, :), 1);
    RawTrace1Hz = RawTrace1Hz(:);
    Bkg1Hz = mean(DataROIFilt1Hz(BkgIdx, :), 1);
    Bkg1Hz = Bkg1Hz(:);

    % Subtract background from the soma trace (1 Hz filtered data)
    Trace1Hz = RawTrace1Hz - Bkg1Hz;
    Trace1Hz = Trace1Hz - mean(Trace1Hz(:));

    % Estimate noise amplitude using negative values in the trace (1 Hz filtered data)
    Mask1Hz = (Trace1Hz < 0);
    NoiseAmp1Hz = sqrt(sum(Trace1Hz.^2 .* Mask1Hz) / sum(Mask1Hz));

    % Compute SNR-normalized raw trace (1 Hz filtered data)
    SNRRawTrace1Hz = Trace1Hz / NoiseAmp1Hz;

    % Store SNRRawTrace1Hz in SpikeInfo
    SpikeInfo.SNRRawTrace1Hz = SNRRawTrace1Hz;

    % -------------------- Rest of the code remains the same --------------------

    % The plotting and additional processing code remains the same.
    % Since plots are generated in the main loop, no further plotting is needed here.
end

%% Function to apply a high-pass filter to video data

function VideoFilt = highpass_video_filt(Param, Video)
    CutOffFreq = Param.CutOffFreq;   % Cutoff frequency for the high-pass filter
    SampleRate = Param.SampleRate;   % Sampling rate of the video data

    NormFreq = CutOffFreq / (SampleRate / 2); % Normalized cutoff frequency
    [bb, aa] = butter(5, NormFreq, 'high');   % Design a 3rd-order Butterworth high-pass filter

    Video = permute(Video, [3 1 2]);          % Rearrange dimensions to apply filter along time
    VideoFilt = filtfilt(bb, aa, Video);      % Apply zero-phase filtering to avoid phase shift
    VideoFilt = permute(VideoFilt, [2 3 1]);  % Restore the original dimensions
end