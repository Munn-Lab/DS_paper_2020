function [SWR, numSWR, freqSWR, fltrdEEG, ripple_locs] = detectSWR_wavelet(eegData,plotfig)

    eeg = eegData{1,1};
    Fs = eegData{1,2};
    time = (1:size(eeg,1))/Fs; % create time vector in seconds
    totalSec = time(end); % total time in seconds
    % a multisite average approach was used to detect SWRs
    swrPassband = [150 250]; % passband range in Hz
%     fltrdEEG = bandpass(eeg, swrPassband(1),swrPassband(2), Fs); %LFPs from all available CA1 cell layer tetrodes were filtered between 150-250Hz
%     
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