function [SWR, numSWR, freqSWR, fltrdEEG, ripple_locs, num_spikes_during_swr,percent_spikes_during_swr] = detectSWR_v2(eegData,plotfig,cellTS)

eeg = eegData{1,1};
Fs = eegData{1,2};
time = (1:size(eeg,1))/Fs; % create time vector in seconds
totalSec = time(end); % total time in seconds
% a multisite average approach was used to detect SWRs
swrPassband = [150 250]; % passband range in Hz
%     fltrdEEG = bandpass(eeg, swrPassband(1),swrPassband(2), Fs); %LFPs from all available CA1 cell layer tetrodes were filtered between 150-250Hz

wOutputs = wFilter(eeg,Fs,200,50,0); %% Changed by Rob to use wavelet-based filtering (200Hz centered wavelet with a bandwidth of 50 (150-250Hz);
fltrdEEG = wOutputs{1};
squaredfltrdEEG = fltrdEEG.^2;% then squared
sigma = 0.004*Fs;
gaussEEG = smoothdata(squaredfltrdEEG,'gaussian', sigma); % This  was smoothed with a Gaussian kernel (sigma=4ms)
analysisEEG = gaussEEG.^(0.5); %the square root of the smoothed sum was analysed

% SWRs were detected when the signal exceeded 3 s.d. of the recording epoch mean for at least 15ms.
stdThreshold = 3;
swrMinLength = 0.015;%15ms
minSWRsampleLength = swrMinLength*Fs;
stddev = std(analysisEEG);
significantRipple = (analysisEEG) > stdThreshold*stddev;  % detected when the signal exceeded 3 s.d.
ripple_locs = significantRipple;
labeledRuns = bwlabel(significantRipple);
measurements = regionprops(labeledRuns, 'Area', 'PixelIdxList');
ripple = analysisEEG(ripple_locs);

for g = 1:length(cellTS) %%% here we're finding the closest value to 0 when we subtract the time of the spike (cellTS) from the sample of EEG (phase_timestamp)
    diff = cellTS(g) - time;
    [~,min_time(g)] = min(abs(diff)); %% This finds the value closest to 0 and returns an index and the value. Here we're supressing the value (by using '~') and getting the index
end

if  ~exist('min_time','var')
spikes_vs_swr = NaN;
cellTS_ripple = NaN;
num_spikes_during_swr = NaN;
percent_spikes_during_swr = NaN;
else
spikes_vs_swr = significantRipple(min_time);
cellTS_ripple = cellTS(spikes_vs_swr);
num_spikes_during_swr = sum(spikes_vs_swr);
percent_spikes_during_swr = num_spikes_during_swr/length(cellTS);
end
SWR = measurements([measurements.Area]>minSWRsampleLength); % at least 15ms.
numSWR = size(SWR,1);
freqSWR = numSWR/totalSec; % number of SWR per Sec
% Plotting SWRs
close all;
% plot random SWR
randomSWR = round(rand(1)*size(SWR,1));

if plotfig == 1
    
    figure(1)
    plot(time(SWR(randomSWR).PixelIdxList), fltrdEEG(SWR(randomSWR).PixelIdxList));
    
    % plot SWR on top of passbanded eeg
    figure(2)
    hold on;
    plot(time, fltrdEEG,'black');
    for i = 1:numSWR
        plot(time(SWR(i).PixelIdxList), fltrdEEG(SWR(i).PixelIdxList),'r');
    end
    
    % plot SWR on top of raw eeg
    figure(3)
    hold on;
    plot(time, eeg,'black');
    for i = 1:numSWR
        plot(time(SWR(i).PixelIdxList), eeg(SWR(i).PixelIdxList),'r');
    end
    
    
    figure(4)
    hold on;
    plot(time, analysisEEG,'black');
    for i = 1:numSWR
        plot(time(SWR(i).PixelIdxList), analysisEEG(SWR(i).PixelIdxList),'r');
    end
else
end
end