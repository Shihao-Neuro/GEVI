% Clear the workspace and close all figures
clear all;
close all;

% Add paths to the subfunctions and example data directories
addpath("./subfunctions/");

% Define the block size for PCA analysis
BNxy = 8; % Block size in x and y directions for PCA analysis

% Define the list of image numbers to process
image_numbers = 24:30; % Modify this to include all the images you wish to process, e.g. 30s1.71camera024-30s1.71camera030

% Loop over each image number
for XX_num = image_numbers
    % Format XX as a zero-padded three-digit string, e.g., '010', '011', etc.
    XX = sprintf('%03d', XX_num);
    
    % Define the file paths and parameters
    FilePath = ['Your_file_path_containing_images_each_in_a_folder\', XX, '\']; % The path where the files are located
    FileNamePrefix = ['30s1.71camera', XX]; % The prefix of the file name, not including the three digits e.g. 001,002 etc
    FileNameSuffix = ['.nd2']; % ND2 file extension
    
    %% Generate motion correction template and estimate motion shifts
    % Read the image stack from an ND2 file using Bio-Formats
    % Ensure that the Bio-Formats toolbox is installed and added to the MATLAB path
    data = bfopen([FilePath, FileNamePrefix, FileNameSuffix]);
    
    % 'data{1}' contains the series data
    series = data{1}; % First series (assuming only one series)
    
    numFrames = size(series,1);
    
    % Get the size of each image
    imageSize = size(series{1,1});
    
    % Preallocate ImStack
    ImStack = zeros(imageSize(1), imageSize(2), numFrames, 'like', series{1,1});
    
    % Stack images into ImStack
    for i = 1:numFrames
        ImStack(:,:,i) = series{i,1};
    end
    
    % Create a template by computing the mean of the image stack across frames
    Template = mean(single(ImStack),3);
    
    %%
    % Define the Region of Interest (ROI) position as [y_start x_start y_size x_size]
    ROIPos = [1 1 size(ImStack,1) size(ImStack,2)]; %[y0 x0 y_width x_width]
    
    % Estimate motion shifts between the template and each frame in the image stack
    Shift = motion_est(Template, ROIPos, ImStack, FilePath, XX);
    
    % Save the template, ROI position, and estimated shifts
    save([FilePath, 'motion_estimation_', XX, '.mat'],'Template','ROIPos','Shift');
    
    %% Apply motion correction
    % Correct the image stack using the estimated shifts
    ImStackCrt = motion_crt(ImStack, Shift);
    
    % Sum adjacent frames to enhance the signal (this can be adjusted as needed)
    ImStackCrt = ImStackCrt(:,:,2:end) + ImStackCrt(:,:,1:end-1);
    
    % Save the corrected image stack
    save([FilePath, '30s1.71camera', XX, '.mat'],'ImStackCrt','-v7.3');
    
    % Prepare for PCA by initializing an empty array
    ImPCA = [];
    
    % Crop the image stack to have dimensions that are multiples of BNxy
    Tmp = ImStackCrt(1:floor(size(ImStackCrt,1)/BNxy)*BNxy,1:floor(size(ImStackCrt,2)/BNxy)*BNxy,:);
    
    % Reshape and average blocks for dimensionality reduction before PCA
    Tmp = reshape(Tmp, BNxy, size(Tmp,1)/BNxy, BNxy, size(Tmp,2)/BNxy, []);
    ImPCA(:,:,:) = squeeze(mean(mean(single(Tmp),1),3));
    
    %%
    % Extract noise components using PCA
    NoisePCA = noise_pca_extract(ImPCA, FilePath, XX);
    
    % Save the PCA results, block size, and processed image data
    save([FilePath, 'pca_', XX, '.mat'],'NoisePCA','BNxy','ImPCA');
    
    % Display progress
    disp(['Processing of image ', XX, ' completed.']);
end

disp('All images processed successfully.');


%% Function to extract noise components using PCA
function NoisePCA = noise_pca_extract(Img, FilePath, XX)
    CmpN = 50; % Number of principal components to compute

    % Rearrange dimensions to prepare for PCA
    Img = permute(Img, [1 2 4 3]);
    Img = reshape(Img, size(Img,1), [], size(Img,4));
    ImgSize = size(Img);
    Img = reshape(Img, ImgSize(1)*ImgSize(2), []);

    %%
    tic
    % Perform Singular Value Decomposition (SVD) for PCA
    [U, S, V] = svds(Img, CmpN);
    toc

    %%
    % Reshape U to get spatial patterns
    Patterns = reshape(U, ImgSize(1), ImgSize(2), []);

    % Compute temporal traces
    Trace = S * V';
    clear S V

    %%
    % Define the number of components to visualize
    numComponentsToPlot = 10; % Adjusted to only plot the first 10 components

    % Create a folder to save the plots
    plotFolder = [FilePath, 'PCA_Plots\'];
    if ~exist(plotFolder, 'dir')
        mkdir(plotFolder);
    end

    % Visualization of the PCA components and saving the plots
    for ii = 1:numComponentsToPlot  % Only loop over the first 10 components
        hFig = figure('Visible', 'off');
        
        % Ensure the figure matches the original appearance
        clf;
        subplot(1,2,1);
        imagesc(Patterns(:,:,ii));
        axis image;
        colorbar;
        title(['Spatial Pattern ' num2str(ii)]);
        
        subplot(1,2,2);
        plot(Trace(ii,:));
        title(['Temporal Trace ' num2str(ii)]);
        
        % Save the figure as an image file
        saveas(hFig, [plotFolder, 'PCA_Component_', num2str(ii), '_', XX, '.png']);
        close(hFig);
    end

    % Store the results in a structure
    NoisePCA.Patterns = Patterns;
    NoisePCA.Trace = Trace;
end


%% Function to estimate motion shifts
function Shift = motion_est(Template, ROIPos, ImStack, FilePath, XX)
    %%
    MaxShift = 6; % Maximum allowed shift in pixels

    % Extract the ROI from the image stack
    ROI = ImStack(ROIPos(1):ROIPos(1)+ROIPos(3)-1, ROIPos(2):ROIPos(2)+ROIPos(4)-1, :);
    ROI = double(ROI);

    % Extract the ROI from the template
    TemplateROI = Template(ROIPos(1):ROIPos(1)+ROIPos(3)-1, ROIPos(2):ROIPos(2)+ROIPos(4)-1, :);

    % Create a mask to ignore the edges where shifts could cause issues
    Mask = ROI(:,:,1)*0;
    Mask(MaxShift+1:end-MaxShift, MaxShift+1:end-MaxShift) = 1;

    % Compute the gradients in the x and y directions
    ROI = (ROI - circshift(ROI, [1 0 0])) + 1i * (ROI - circshift(ROI, [0 1 0]));
    TemplateROI = (TemplateROI - circshift(TemplateROI, [1 0])) + 1i * (TemplateROI - circshift(TemplateROI, [0 1]));

    % Apply the mask to the ROI and template ROI
    ROIRaw = ROI .* Mask;
    TemplateROI = TemplateROI .* Mask;

    % Compute the cross-correlation between the template ROI and each frame's ROI
    XCorr = abs(fftshift(fftshift(ifft2(fft2(TemplateROI) .* conj(fft2(ROIRaw))),1),2));

    % Compute normalization factors for the cross-correlation
    Ref = sqrt(abs(fftshift(fftshift(ifft2(fft2(Mask) .* conj(fft2(abs(ROIRaw).^2))),1),2)));
    Ref = Ref .* sqrt(abs(fftshift(fftshift(ifft2(fft2(abs(TemplateROI).^2) .* conj(fft2(Mask))),1),2)));

    % Normalize the cross-correlation to obtain correlation coefficients
    XCorr = XCorr ./ (Ref + eps);

    % Determine the center of the cross-correlation
    Center = floor(size(XCorr) / 2 + 1);

    % Extract a region around the center corresponding to feasible shifts
    XCorrROI = XCorr(Center(1)-MaxShift:Center(1)+MaxShift, Center(2)-MaxShift:Center(2)+MaxShift, :);
    Tmp = reshape(XCorrROI, [], size(XCorrROI,3));

    % Find the peak correlation value and its index for each frame
    [XCorrMax, Idx] = max(Tmp, [], 1);
    [YIdx, XIdx] = ind2sub([size(XCorrROI,1), size(XCorrROI,2)], Idx);

    % **Calculate the coarse shifts relative to the center position**
    YShift = YIdx - MaxShift - 1;
    XShift = XIdx - MaxShift - 1;

    % Define the size of the area around the peak for sub-pixel shift estimation
    PeakFitSize = 2;

    % Initialize arrays for sub-pixel shifts
    YIdxFine = zeros(1, size(XCorr,3));
    XIdxFine = zeros(1, size(XCorr,3));

    % Compute the sub-pixel shifts by fitting a centroid around the peak
    for ii = 1:size(XCorr,3)
        % Calculate the valid range for Y and X indices
        y_start = max(YIdx(ii) - PeakFitSize, 1);
        y_end = min(YIdx(ii) + PeakFitSize, size(XCorrROI,1));
        x_start = max(XIdx(ii) - PeakFitSize, 1);
        x_end = min(XIdx(ii) + PeakFitSize, size(XCorrROI,2));

        % Adjust the meshgrid accordingly
        [xx_tmp, yy_tmp] = meshgrid(x_start - XIdx(ii):x_end - XIdx(ii), y_start - YIdx(ii):y_end - YIdx(ii));

        % Extract the neighborhood around the peak
        Tmp = XCorrROI(y_start:y_end, x_start:x_end, ii);
        Tmp = Tmp - min(Tmp(:));

        % Compute sub-pixel shifts using centroid method
        YIdxFine(ii) = sum(sum(yy_tmp .* Tmp)) / sum(Tmp(:)) + YIdx(ii);
        XIdxFine(ii) = sum(sum(xx_tmp .* Tmp)) / sum(Tmp(:)) + XIdx(ii);
    end

    % Calculate the shifts relative to the center position
    YShiftFine = YIdxFine - MaxShift - 1;
    XShiftFine = XIdxFine - MaxShift - 1;

    % Prepare the data for plotting
    YShift = YShift(:); % Ensure column vector
    XShift = XShift(:);
    YShiftFine = YShiftFine(:);
    XShiftFine = XShiftFine(:);

    % Check that all shift vectors are the same length
    numFrames = length(YShift);
    assert(length(XShift) == numFrames && length(YShiftFine) == numFrames && length(XShiftFine) == numFrames, 'Shift vectors must be the same length.');

    % Optionally plotting of the estimated shifts and correlation maxima
    % Plot both coarse and fine shifts to match the original code
    hFigShifts = figure('Visible', 'off');

    % Create a matrix with each shift as a column
    shiftMatrix = [YShift, XShift, YShiftFine, XShiftFine];

    % Plot the shift data
    plot(shiftMatrix);

    % Add labels and legend
    legend('Y Shift', 'X Shift', 'Y Shift (Fine)', 'X Shift (Fine)');
    title('Estimated Motion Shifts');
    ylabel('Shift (pixels)');
    xlabel('Frame Number');

    % Save and close the figure
    saveas(hFigShifts, [FilePath, 'Estimated_Motion_Shifts_', XX, '.png']);
    close(hFigShifts);

    % Plot the maximum cross-correlation values
    hFigXCorr = figure('Visible', 'off');
    plot(XCorrMax);
    title('Maximum Cross-Correlation Values');
    ylabel('Cross-Correlation');
    xlabel('Frame Number');
    saveas(hFigXCorr, [FilePath, 'XCorr_Maxima_', XX, '.png']);
    close(hFigXCorr);

    % Store the fine shifts in a structure
    Shift.YShift = YShiftFine;
    Shift.XShift = XShiftFine;

    % Optionally, you may also store the coarse shifts if needed
    Shift.YShiftCoarse = YShift;
    Shift.XShiftCoarse = XShift;
end


%%
% ******************Function to correct the images using the estimated shifts
function ImgCrt = motion_crt(Img, Shift)
    MaxShift = 6; % Maximum allowed shift in pixels

    XShiftFine = Shift.XShift;
    YShiftFine = Shift.YShift;
    TN = length(XShiftFine); % Total number of frames to correct

    % Initialize the corrected image stack
    ImgCrt = double(Img(:,:,end-TN+1:end));

    % Create coordinate grids for the images
    [x0, y0] = meshgrid(1:size(Img,2), 1:size(Img,1));

    % Apply the calculated shifts to each frame using interpolation
    for ii = 1:length(XShiftFine)
        ImgCrt(:,:,ii) = interp2(x0, y0, ImgCrt(:,:,ii), x0 - XShiftFine(ii), y0 - YShiftFine(ii), 'linear', 0);
    end

    % Create a mask to exclude regions affected by the shifts
    Mask = ImgCrt(:,:,1) * 0;
    Mask(MaxShift+1:end-MaxShift, MaxShift+1:end-MaxShift) = 1;

    % Apply the mask to the corrected images and convert to uint16 format
    ImgCrt = uint16(ImgCrt .* Mask);
end
