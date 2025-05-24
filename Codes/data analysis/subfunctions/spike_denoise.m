function SpikeInfo = spike_denoise(Param, RawTrace)
    SampleRate = Param.SampleRate;
    SpikeTemplateLength = Param.SpikeTemplateLength;
    SpikeTemplateN = Param.SpikeTemplateN;
    CutOffFreq = Param.CutOffFreq;
    HeadTailSize = Param.HeadTailSize;
    SNRList = Param.SNRList;
    MinSpikeTemplateN = Param.MinSpikeTemplateN;

    %%
    DataHP = highpass_filt(RawTrace, CutOffFreq, SampleRate);
    % Capture SNRThdTheory from template_threshold function
    [SpikeTemplateIdx, SNRThd, PosNegPeakCnt, SNRThdTheory] = template_threshold(DataHP, HeadTailSize, SNRList, MinSpikeTemplateN);

    if ~isempty(SpikeTemplateIdx)
        [SpikeIdx, RawTraceNoiseAmp] = simple_threshold(DataHP, SNRThd);

        SpikeTemplate = spike_template_gen(DataHP, SpikeTemplateIdx, SpikeTemplateLength, SpikeTemplateN);
        DataFilt = whitened_matched_filter(DataHP, SpikeTemplateIdx, SpikeTemplateLength, SpikeTemplateN);
        [SpikeIdx, FiltTraceNoiseAmp] = simple_threshold(DataFilt, SNRThd);
        SpikeRecon = DataFilt * 0;
        SpikeRecon(SpikeIdx) = 1;
        SpikeRecon = conv(SpikeRecon, SpikeTemplate, 'same');
        Coef = sum(SpikeRecon .* RawTrace) / sum(SpikeRecon .^ 2);
        SpikeSub = RawTrace - Coef * SpikeRecon;

        % Store spike information
        SpikeInfo.RawTrace = RawTrace;
        SpikeInfo.FiltTrace = DataFilt;

        SpikeInfo.SpikeIdx = SpikeIdx;
        SpikeInfo.SpikeTemplate = SpikeTemplate;
        SpikeInfo.SpikeRecon = SpikeRecon;
        SpikeInfo.SpikeSub = SpikeSub;
        SpikeInfo.SpikeTemplateIdx = SpikeTemplateIdx;
        SpikeInfo.SNRThd = SNRThd;
        SpikeInfo.PosNegPeakCnt = PosNegPeakCnt;
        % Store SNRThdTheory
        SpikeInfo.SNRThdTheory = SNRThdTheory;

        % Calculate and store all spike amplitudes
        SpikeRawAmp = RawTrace(SpikeIdx);
        SpikeFiltAmp = DataFilt(SpikeIdx);
        % Store unsorted amplitudes
        SpikeInfo.SpikeRawAmp = SpikeRawAmp;
        SpikeInfo.SpikeFiltAmp = SpikeFiltAmp;

        % For SNR calculation, sort the amplitudes in descending order
        SpikeRawAmpSorted = sort(SpikeRawAmp, 'descend');
        SpikeFiltAmpSorted = sort(SpikeFiltAmp, 'descend');

        % Calculate Noise Amplitudes
        MaskRaw = RawTrace < 0;
        NoiseAmpRaw = sqrt(sum(RawTrace .^ 2 .* MaskRaw) / sum(MaskRaw));
        SNRRaw = mean(SpikeRawAmpSorted(1 : min(10, length(SpikeRawAmpSorted)))) / NoiseAmpRaw;

        MaskFilt = DataFilt < 0;
        NoiseAmpFilt = sqrt(sum(DataFilt .^ 2 .* MaskFilt) / sum(MaskFilt));
        SNRFilt = mean(SpikeFiltAmpSorted(1 : min(10, length(SpikeFiltAmpSorted)))) / NoiseAmpFilt;

        SpikeN = length(SpikeIdx);

        % Add additional information to SpikeInfo
        SpikeInfo.SNRRaw = SNRRaw;
        SpikeInfo.SNRFilt = SNRFilt;
        SpikeInfo.NoiseAmpRaw = NoiseAmpRaw;
        SpikeInfo.NoiseAmpFilt = NoiseAmpFilt;
        SpikeInfo.SpikeN = SpikeN;

        % Optional: If you wish to store the sorted amplitudes as well
        SpikeInfo.SpikeRawAmpSorted = SpikeRawAmpSorted;
        SpikeInfo.SpikeFiltAmpSorted = SpikeFiltAmpSorted;

        % Add the SNRmedianRaw, SNRpeakRaw, and SNRmeanRaw to SpikeInfo
        SNRmedianRaw = median(SpikeRawAmpSorted / NoiseAmpRaw);
        SNRpeakRaw = max(SpikeRawAmpSorted / NoiseAmpRaw);
        SNRmeanRaw = mean(SpikeRawAmpSorted / NoiseAmpRaw);
        SpikeInfo.SNRmedianRaw = SNRmedianRaw;
        SpikeInfo.SNRpeakRaw = SNRpeakRaw;
        SpikeInfo.SNRmeanRaw = SNRmeanRaw;

        % Add the SNRmedianFilt, SNRpeakFilt, and SNRmeanFilt to SpikeInfo
        SNRmedianFilt = median(SpikeFiltAmpSorted / NoiseAmpFilt);
        SNRpeakFilt = max(SpikeFiltAmpSorted / NoiseAmpFilt);
        SNRmeanFilt = mean(SpikeFiltAmpSorted / NoiseAmpFilt);
        SpikeInfo.SNRmedianFilt = SNRmedianFilt;
        SpikeInfo.SNRpeakFilt = SNRpeakFilt;
        SpikeInfo.SNRmeanFilt = SNRmeanFilt;

        % Compute the SNR-normalized traces
        SpikeInfo.SNRRawTrace = RawTrace / NoiseAmpRaw;
        SpikeInfo.SNRFiltTrace = DataFilt / NoiseAmpFilt;

        % Compute SNRRawTrace1Hz
        DataHP_1Hz = highpass_filt(RawTrace, 1, SampleRate);

        % Compute NoiseAmpRaw1Hz
        MaskRaw1Hz = DataHP_1Hz < 0;
        NoiseAmpRaw1Hz = sqrt(sum(DataHP_1Hz .^ 2 .* MaskRaw1Hz) / sum(MaskRaw1Hz));
        SpikeInfo.NoiseAmpRaw1Hz = NoiseAmpRaw1Hz; % Store in SpikeInfo

        % Compute SNRRawTrace1Hz
        SNRRawTrace1Hz = DataHP_1Hz / NoiseAmpRaw1Hz;

        % Store in SpikeInfo
        SpikeInfo.SNRRawTrace1Hz = SNRRawTrace1Hz;

        % Compute SNRmedianRaw1Hz, SNRpeakRaw1Hz, SNRmeanRaw1Hz
        SpikeAmpRaw1Hz = DataHP_1Hz(SpikeIdx);
        SNRvaluesRaw1Hz = SpikeAmpRaw1Hz / NoiseAmpRaw;  % Note: Using NoiseAmpRaw as per your request
        SNRmedianRaw1Hz = median(SNRvaluesRaw1Hz);
        SNRpeakRaw1Hz = max(SNRvaluesRaw1Hz);
        SNRmeanRaw1Hz = mean(SNRvaluesRaw1Hz);

        % Store these values in SpikeInfo
        SpikeInfo.SNRmedianRaw1Hz = SNRmedianRaw1Hz;
        SpikeInfo.SNRpeakRaw1Hz = SNRpeakRaw1Hz;
        SpikeInfo.SNRmeanRaw1Hz = SNRmeanRaw1Hz;

        % Additions to compute and store mean ± SEM spike waveforms

        % Define parameters for spike waveform extraction
        pre_points = 5; % Points before the spike index
        post_points = 5; % Points after the spike index

        % Extract spike waveforms from SNRRawTrace1Hz
        SpikeWaveforms1Hz = [];
        for i = 1:length(SpikeIdx)
            idx = SpikeIdx(i);
            % Check if indices are within bounds
            if (idx - pre_points >= 1) && (idx + post_points <= length(SNRRawTrace1Hz))
                snippet = SNRRawTrace1Hz((idx - pre_points):(idx + post_points));
                SpikeWaveforms1Hz = [SpikeWaveforms1Hz; snippet'];
            end
        end

        % Compute mean and SEM for SNRRawTrace1Hz
        if ~isempty(SpikeWaveforms1Hz)
            mean_waveform_1Hz = mean(SpikeWaveforms1Hz, 1);
            sem_waveform_1Hz = std(SpikeWaveforms1Hz, 0, 1) / sqrt(size(SpikeWaveforms1Hz, 1));
        else
            mean_waveform_1Hz = [];
            sem_waveform_1Hz = [];
        end

        % Store in SpikeInfo
        SpikeInfo.Waveform_Raw1Hz.mean = mean_waveform_1Hz;
        SpikeInfo.Waveform_Raw1Hz.sem = sem_waveform_1Hz;

        % Extract spike waveforms from SNRRawTrace
        SpikeWaveformsRaw = [];
        for i = 1:length(SpikeIdx)
            idx = SpikeIdx(i);
            % Check if indices are within bounds
            if (idx - pre_points >= 1) && (idx + post_points <= length(SpikeInfo.SNRRawTrace))
                snippet = SpikeInfo.SNRRawTrace((idx - pre_points):(idx + post_points));
                SpikeWaveformsRaw = [SpikeWaveformsRaw; snippet'];
            end
        end

        % Compute mean and SEM for SNRRawTrace
        if ~isempty(SpikeWaveformsRaw)
            mean_waveform_Raw = mean(SpikeWaveformsRaw, 1);
            sem_waveform_Raw = std(SpikeWaveformsRaw, 0, 1) / sqrt(size(SpikeWaveformsRaw, 1));
        else
            mean_waveform_Raw = [];
            sem_waveform_Raw = [];
        end

        % Store in SpikeInfo
        SpikeInfo.Waveform_Raw25Hz.mean = mean_waveform_Raw;
        SpikeInfo.Waveform_Raw25Hz.sem = sem_waveform_Raw;

        % Extract spike waveforms from SNRFiltTrace
        SpikeWaveformsFilt = [];
        for i = 1:length(SpikeIdx)
            idx = SpikeIdx(i);
            % Check if indices are within bounds
            if (idx - pre_points >= 1) && (idx + post_points <= length(SpikeInfo.SNRFiltTrace))
                snippet = SpikeInfo.SNRFiltTrace((idx - pre_points):(idx + post_points));
                SpikeWaveformsFilt = [SpikeWaveformsFilt; snippet'];
            end
        end

        % Compute mean and SEM for SNRFiltTrace
        if ~isempty(SpikeWaveformsFilt)
            mean_waveform_Filt = mean(SpikeWaveformsFilt, 1);
            sem_waveform_Filt = std(SpikeWaveformsFilt, 0, 1) / sqrt(size(SpikeWaveformsFilt, 1));
        else
            mean_waveform_Filt = [];
            sem_waveform_Filt = [];
        end

        % Store in SpikeInfo
        SpikeInfo.Waveform_Filt.mean = mean_waveform_Filt;
        SpikeInfo.Waveform_Filt.sem = sem_waveform_Filt;

        % Optional: Plotting the mean ± SEM waveforms

        % Plot for SNRRawTrace1Hz
        if ~isempty(mean_waveform_1Hz)
            figure(1002);
            t = (-pre_points):post_points; % Time axis relative to spike peak
            % Plot mean waveform
            plot(t, mean_waveform_1Hz, 'b', 'LineWidth', 2);
            hold on;
            % Plot mean ± SEM area
            fill([t, fliplr(t)], [mean_waveform_1Hz + sem_waveform_1Hz, fliplr(mean_waveform_1Hz - sem_waveform_1Hz)], ...
                'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            xlabel('Time (frames)');
            ylabel('Amplitude');
            title('Mean ± SEM Spike Waveform (SNRRawTrace1Hz)');
            hold off;
        else
            warning('No valid spike waveforms found for plotting (SNRRawTrace1Hz).');
        end

        % Plot for SNRRawTrace
        if ~isempty(mean_waveform_Raw)
            figure(1003);
            t = (-pre_points):post_points; % Time axis relative to spike peak
            % Plot mean waveform
            plot(t, mean_waveform_Raw, 'g', 'LineWidth', 2);
            hold on;
            % Plot mean ± SEM area
            fill([t, fliplr(t)], [mean_waveform_Raw + sem_waveform_Raw, fliplr(mean_waveform_Raw - sem_waveform_Raw)], ...
                'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            xlabel('Time (frames)');
            ylabel('Amplitude');
            title('Mean ± SEM Spike Waveform (SNRRawTrace)');
            hold off;
        else
            warning('No valid spike waveforms found for plotting (SNRRawTrace).');
        end

        % Plot for SNRFiltTrace
        if ~isempty(mean_waveform_Filt)
            figure(1004);
            t = (-pre_points):post_points; % Time axis relative to spike peak
            % Plot mean waveform
            plot(t, mean_waveform_Filt, 'r', 'LineWidth', 2);
            hold on;
            % Plot mean ± SEM area
            fill([t, fliplr(t)], [mean_waveform_Filt + sem_waveform_Filt, fliplr(mean_waveform_Filt - sem_waveform_Filt)], ...
                'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            xlabel('Time (frames)');
            ylabel('Amplitude');
            title('Mean ± SEM Spike Waveform (SNRFiltTrace)');
            hold off;
        else
            warning('No valid spike waveforms found for plotting (SNRFiltTrace).');
        end

        % Plotting raw, filtered, and reconstructed traces
        figure(2000);
        ax1 = subplot(3,1,1); plot(RawTrace); title('Raw Trace');
        ax2 = subplot(3,1,2); plot(DataFilt); title('Filtered Trace');
        ax3 = subplot(3,1,3); plot(SpikeRecon); title('Spike Reconstruction');
        linkaxes([ax1, ax2, ax3], 'x');
    else
        DataFilt = whitened_matched_filter(DataHP, SpikeTemplateIdx, SpikeTemplateLength, SpikeTemplateN);

        % Compute Noise Amplitudes
        MaskRaw = RawTrace < 0;
        NoiseAmpRaw = sqrt(sum(RawTrace .^ 2 .* MaskRaw) / sum(MaskRaw));

        MaskFilt = DataFilt < 0;
        NoiseAmpFilt = sqrt(sum(DataFilt .^ 2 .* MaskFilt) / sum(MaskFilt));

        % Compute the SNR-normalized traces
        SpikeInfo.SNRRawTrace = RawTrace / NoiseAmpRaw;
        SpikeInfo.SNRFiltTrace = DataFilt / NoiseAmpFilt;

        % Compute SNRRawTrace1Hz
        DataHP_1Hz = highpass_filt(RawTrace, 1, SampleRate);

        MaskRaw1Hz = DataHP_1Hz < 0;
        NoiseAmpRaw1Hz = sqrt(sum(DataHP_1Hz .^ 2 .* MaskRaw1Hz) / sum(MaskRaw1Hz));
        SpikeInfo.NoiseAmpRaw1Hz = NoiseAmpRaw1Hz; % Store in SpikeInfo

        SNRRawTrace1Hz = DataHP_1Hz / NoiseAmpRaw1Hz;

        % Store in SpikeInfo
        SpikeInfo.SNRRawTrace1Hz = SNRRawTrace1Hz;

        % Store SNRThdTheory in SpikeInfo
        SpikeInfo.SNRThdTheory = SNRThdTheory;

        % Set other fields of SpikeInfo to empty
        SpikeInfo.RawTrace = RawTrace;
        SpikeInfo.FiltTrace = DataFilt;

        SpikeInfo.SpikeIdx = [];
        SpikeInfo.SpikeTemplate = [];
        SpikeInfo.SpikeRecon = [];
        SpikeInfo.SpikeSub = [];
        SpikeInfo.SpikeTemplateIdx = SpikeTemplateIdx;
        SpikeInfo.SNRThd = SNRThd;
        SpikeInfo.PosNegPeakCnt = PosNegPeakCnt;

        SpikeInfo.SpikeRawAmp = [];
        SpikeInfo.SpikeFiltAmp = [];
        SpikeInfo.SpikeRawAmpSorted = [];
        SpikeInfo.SpikeFiltAmpSorted = [];
        SpikeInfo.SNRRaw = [];
        SpikeInfo.SNRFilt = [];
        SpikeInfo.NoiseAmpRaw = NoiseAmpRaw;
        SpikeInfo.NoiseAmpFilt = NoiseAmpFilt;
        SpikeInfo.SpikeN = 0;

        % Set SNRmedianRaw, SNRpeakRaw, and SNRmeanRaw to empty
        SpikeInfo.SNRmedianRaw = [];
        SpikeInfo.SNRpeakRaw = [];
        SpikeInfo.SNRmeanRaw = [];

        % Set SNRmedianFilt, SNRpeakFilt, and SNRmeanFilt to empty
        SpikeInfo.SNRmedianFilt = [];
        SpikeInfo.SNRpeakFilt = [];
        SpikeInfo.SNRmeanFilt = [];

        % Set SNRmedianRaw1Hz, SNRpeakRaw1Hz, SNRmeanRaw1Hz to empty
        SpikeInfo.SNRmedianRaw1Hz = [];
        SpikeInfo.SNRpeakRaw1Hz = [];
        SpikeInfo.SNRmeanRaw1Hz = [];

        % Set spike waveform fields to empty
        SpikeInfo.Waveform_Raw1Hz.mean = [];
        SpikeInfo.Waveform_Raw1Hz.sem = [];
        SpikeInfo.Waveform_Raw25Hz.mean = [];
        SpikeInfo.Waveform_Raw25Hz.sem = [];
        SpikeInfo.Waveform_Filt.mean = [];
        SpikeInfo.Waveform_Filt.sem = [];

        % Plotting raw and filtered traces
        figure(2000);
        ax1 = subplot(3,1,1); plot(RawTrace); title('Raw Trace');
        ax2 = subplot(3,1,2); plot(DataFilt); title('Filtered Trace');
        ax3 = subplot(3,1,3); plot(zeros(size(RawTrace))); title('Spike Reconstruction');
        linkaxes([ax1, ax2, ax3], 'x');
    end
end



%%
function DataOut = whitened_matched_filter(DataIn,SpikeIdx, SpikeTemplateLength,SpikeTemplateN)
    SpikeRecon=DataIn*0;
    SpikeRecon(SpikeIdx)=1;
    SpikeMask=conv(SpikeRecon,ones(2*SpikeTemplateLength+1,1),'same');
    NoiseMask=SpikeMask<0.5;
    
    Noise=DataIn(find(NoiseMask));
    NoisePSD=pwelch(Noise,1000,[],2*length(DataIn)-1);
    NoisePSD=[NoisePSD(:);flipud(NoisePSD(1:end-1))];
    DataWhiten=real(ifft(fft(DataIn,length(NoisePSD))./sqrt(NoisePSD)));
    DataWhiten=DataWhiten(1:length(DataIn));
    SpikeTemplate=spike_template_gen(DataWhiten,SpikeIdx,SpikeTemplateLength,SpikeTemplateN);
    DataOut=conv(DataWhiten,flipud(SpikeTemplate),'same');
end

function SpikeTemplate=spike_template_gen(Data,SpikeIdx,SpikeTemplateLength,SpikeTemplateN)
    SpikeIdxMask=(SpikeIdx>SpikeTemplateLength).*(SpikeIdx<(length(Data)-SpikeTemplateLength));
    SpikeIdx=SpikeIdx(find(SpikeIdxMask));
    
    if length(SpikeIdx)>SpikeTemplateN
        [~,TmpIdx]=sort(Data(SpikeIdx),'descend');
        SpikeIdx=SpikeIdx(TmpIdx(1:SpikeTemplateN));
    end
    
    SpikeTemplate=[];
    for ii=1:length(SpikeIdx)
        SpikeTemplate(:,ii)=Data(SpikeIdx(ii)-SpikeTemplateLength:SpikeIdx(ii)+SpikeTemplateLength);
    end
    
    % figure(100);
    % plot(SpikeTemplate,'color', [0.5 0.5 0.5]);hold on;
    SpikeTemplate=mean(SpikeTemplate,2);
    % plot(SpikeTemplate,'k','linewidth',2);hold off;
end

function [PeakIdx, NoiseAmp]=simple_threshold(Data,SNRThd)
    Mask=Data<0;
    NoiseAmp=sqrt(sum(Mask.*Data.^2)/sum(Mask));
    Data=Data.*(Data>(SNRThd*NoiseAmp));
    [~,PeakIdx]=findpeaks(Data);
end

%%
function [PeakTemplateIdx, SNRThd, PosNegPeakCnt, SNRThdTheory] = template_threshold(Data, HeadTailSize, SNRList, MinSpikeTemplateN)
    %% Spike must be a single peak with 2 preceding and following data points smaller than the peak
    
    Data = Data(HeadTailSize + 1:end - HeadTailSize);
    Mask = Data < 0;
    NoiseAmp = sqrt(sum(Mask .* Data .^ 2) / sum(Mask));
    Data = Data / NoiseAmp;
    PosPeaks = (Data > circshift(Data, 1)) & (Data > circshift(Data, -1));
    NegPeaks = (Data < circshift(Data, 1)) & (Data < circshift(Data, -1));
    PosPeakIdx = find(PosPeaks);
    NegPeakIdx = find(NegPeaks);
    
    % Compute SNRThdTheory
    SNRThdTheory = 0.5 * (1 + erf(-SNRList / sqrt(2))) * length(Data(:));
    
    SinglePeaks = (Data > circshift(Data, 1)) & (Data > circshift(Data, 2)) & (Data > circshift(Data, -1)) & (Data > circshift(Data, -2));
    SinglePeakIdx = find(SinglePeaks);
    
    PosNegPeakCnt = [];
    PeakTemplateIdx = [];
    for ii = 1:length(SNRList)
        PosNegPeakCnt(ii, 1) = sum(Data(PosPeakIdx) > SNRList(ii));
        PosNegPeakCnt(ii, 2) = sum(Data(NegPeakIdx) < -SNRList(ii));
    end
    
    % Print PosNegPeakCnt
    disp('PosNegPeakCnt:');
    disp(PosNegPeakCnt);
    
    % Conditions to select SNR threshold
    %Tmp1 = PosNegPeakCnt(:, 1) > (max(PosNegPeakCnt(:, 2), SNRThdTheory(:))*3);
    Tmp1 = PosNegPeakCnt(:, 1) > (max(PosNegPeakCnt(:, 2), SNRThdTheory(:)));
    Tmp2 = (PosNegPeakCnt(:, 1) - PosNegPeakCnt(:, 2)) > MinSpikeTemplateN;
    SNRThd = SNRList(find(Tmp1 & Tmp2, 1));
    
    if ~isempty(SNRThd)
        SinglePeaks = (Data > circshift(Data, 1)) & (Data > circshift(Data, 2)) & (Data > circshift(Data, -1)) & (Data > circshift(Data, -2));
        SinglePeakIdx = find(SinglePeaks);
        TmpIdx = find(Data(SinglePeakIdx) > SNRThd);
        PeakTemplateIdx = SinglePeakIdx(TmpIdx) + HeadTailSize;
    end
end


function DataOut=highpass_filt(DataIn,CutOffFreq,SampleRate)
    NormFreq=CutOffFreq/(SampleRate/2);
    [bb,aa]=butter(2,NormFreq,'high');
    DataOut=filtfilt(bb,aa,double(DataIn));
    DataOut=DataOut-median(DataOut(:));
end
