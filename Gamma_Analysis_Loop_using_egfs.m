
%__________________________________________________________________________
%
%                       Program parameters
%__________________________________________________________________________

p_struct
% Name for the folder where the images will be stored. In addition the name
% of the input file will be used in the folder name. Example: in121314.txt
% will give a folder name gammaImages_in121314
p.imageFolder = 'gammaImages';
% Name for folder where result data from the program will be stored.
% In addition the name of the input file will be used in the folder name. 
% Example: in121314.txt will give a folder name gammaResultsData_in121314.
p.dataFolder = 'gammaResultsData';
% Bin width for the rate maps
p.binWidth = 2.5; % [cm]
% Minimum time bins in the rate map must have been visited by the rat. Bins
% with less time will be set to NaN, and plotted as white pixels. Minimum
% possible value is 0.020 for Axona data and 0.040 for NeuraLynx data
% because of the position sampling rate.
p.minBinTime = 0.020; % [sec
% Format for the images
% format = 1 -> bmp (24 bit)
% format = 2 -> png
% format = 3 -> eps
% format = 4 -> jpg
% format = 5 -> ai (Adobe Illustrator)
% format = 6 -> tiff (24 bit)
% format = 7 -> fig (Matlab figure)
% format = 8 -> svg (Vector graphics)
p.imageFormat = 8;
% Frequency vector used in generating some of the figures. Normally it
% should be 2:2:200, which means it goes from 2 to 200 Hz with a resolution
% of 2 Hz.
p.frequencyVector = 2:2:200;
% Width of the Morlet wavelet used in this program. A value of 7 Hz has been
% found to be good. You should never go below 5 Hz, because it may lead to
% unwanted effects at low frequencies.
p.waveletWidth = 7; % [Hz]
% Theta border frequencies that define the theta band (Laura Colgin used 6
% to 10 Hz)
p.thetaFrequency = [6, 10]; 
% Frequency borders for the slow and fast gamma bands. Use default values
% (Slow = [25, 55] Fast = [65, 140]) if you don't now the exact values for
% your EEG.
p.slowGammaFrequency = [25, 55];
p.fastGammaFrequency = [65, 140];
% This program has been optimized for running on multicore CPU's. Different
% parts of the program will make use of parallel execution of calculations. 
% In this parameter you can set the number of cpu cores you have available.
% Setting it to 'max' will make Matlab use the maximum possible for your 
% computer. 
%       Example:    p.numCpuCores = 'max';
%
% If you don't want it to use the maximum you can set it to a lower value
% by typing the number of cpu cores you want it to run on.
%       Example:    p.numCpuCores = 2;
%
% On a local computer the maximum number of cores you can set is 8 even if 
% your computer has more cores, this is because of limitations in the 
% Matlab parallel computing toolbox. If you don't have the parallel
% computing toolbox the code will be run on only one cpu.
p.numCpuCores = 'max';
% Set what analysis to do. A zero mark that analysis will not be done while
% a one mark that it have to be done
p.analysisList = zeros(11,1);
% 1: Cross-frequency coherence
p.analysisList(1) = 1;
% 2: Theta phase of the gamma episodes
p.analysisList(2) = 1;
% Time-frequency representation of the EEG
p.analysisList(3) = 1;
% Make figures of slows and fast gamma windows
p.analysisList(4) = 1;
% Save slow and fast gamma power to file. Can be used to make a figure of
% mean slow and fast gamma power ratios (Figure 2.a in paper). The data
% will be stored in the datafolder as a mat-file. Data will be stored in
% the file sessionId_GammaPower.mat.
p.analysisList(5) = 1;
% Make figure with interregional synchrony during slow and fast gamma
% oscillations.
p.analysisList(6) = 0;
% Make figure of cell firing during fast and slow gamma oscillations and
% save data to file. Data is saved in the file
% sessionId_cellId_GammaPower.mat
p.analysisList(7) = 1;
% Make figures of consecutive theta cycles with slow and fast gamma and
% save data to files.
p.analysisList(8) = 1;
% Set if the rate maps for cell firing during slow and fast gamma have to
% be made.
p.analysisList(9) = 1;
% Make figure of histogram for intervals between slow and fast gamma peaks 
% and cell spike times.
p.analysisList(10) = 1;
% Set if the correlation of slow gamma versus fast gamma have to be done.
% Result will be save in the file sessionId_slowFastGammaCorrelation.mat
p.analysisList(11) = 1;

EGF = table;
[file,path] = uigetfile('*.egf','Choose the .egf files','multiselect','on');
cd(path)
p.imageFormat = 8;
Gamma_Results = table;
if ischar(file)
    for h = 1:size(file,1)
    EGF.path(h,1) = strcat(string(path),file);
    end
else
for h = 1:length(file)
    EGF.path(h,1) = strcat(string(path),file{h});
end
end
plotfig = 1;
for z = 1:height(EGF)
        [status,eeg,samplerate,bytesPerSample] = readEGF(EGF.path(z,1));
        File = strcat(EGF.path(z,1));
        eegData1{1} = eeg; eegData1{2} = samplerate;
        maxTime = (length(eeg)-1) / samplerate;
        eegTS = 0:1/samplerate:maxTime;
        eegData1{4} = eegTS;
        eegData1{5} = 1;
        [CRF,freqVec1,freqVec2] = crossfreqCoh(eegData1{1},eegData1{1},p.frequencyVector ,eegData1{2},p);
        thetaCRF = mean(CRF(:,freqVec1 >= p.thetaFrequency(1) & freqVec1 <= p.thetaFrequency(2)),2);    
        meanCoherence = mean(thetaCRF);
        Gamma_Results.Mean_theta_gamma_coherence(z,1) = meanCoherence;
         % slow gamma
        [slowGammaPeakInd, slowGammaPeakIndNonOverlap, slowGammaStart, slowGammaStop, slowGammaWindowsEEG,] = gammaEpisodes_one_channel(eegData1{1}...
             , p.slowGammaFrequency(1), p.slowGammaFrequency(2), eegData1{2});
         % fast gamma
          [fastGammaPeakInd, fastGammaPeakIndNonOverlap, fastGammaStart, fastGammaStop, fastGammaWindowsEEG,] = gammaEpisodes_one_channel(eegData1{1}...
        , p.fastGammaFrequency(1), p.fastGammaFrequency(2), eegData1{2});
         % find the theta phase of fast v slow theta epochs
          [thetaPhaseBinSlowGamma, thetaPhaseBinFastGamma] = thetaPhaseOfGammaEpisodes(eegData1{1}, slowGammaPeakIndNonOverlap, fastGammaPeakIndNonOverlap, eegData1{2}, p);
        Gamma_Results.Theta_phase_slow_gamma{z} = thetaPhaseBinSlowGamma;
        Gamma_Results.Theta_phase_fast_gamma{z} = thetaPhaseBinFastGamma;
        disp('Calculating the time-frequency representation of the EEG using the Morlet wavelet')
        disp('This takes a few minutes, please wait ...')
        % Calculate the time frequency representation of the eeg (with
        % prewhitening of the signal)
        [thetaWindows, TFR, thetaWindowStart, thetaWindowStop] = TFR_Analysis(eegData1{1}, p.frequencyVector, eegData1{2}, p.waveletWidth);
        
        % Number of theta windows
        % numThetaWindows = size(thetaWindows,2);
    
        % Calculate the mean theta
        meanTheta = mean(thetaWindows,2);
        
        % Locate the trough
        [~, ind] = min(meanTheta);
        
        % Make the time vector for plotting
        timeVector = 1:length(meanTheta);
        timeVector = ((timeVector - ind(1)) / eegData1{2}) * 1000;

        
        numFrequencies = length(p.frequencyVector);
        windowLength = length(timeVector);
        
        % Set the start and stop for the plots
        timeStart = find(timeVector >= -100, 1);
        timeStop = find(timeVector >= 20, 1);
        
        % Set the stop frequency for the plots
        freqStop = find(p.frequencyVector >= 140, 1);
        
        % Calculate the mean TFR over all theta cycles
        meanTFR = mean(TFR,1);
        meanTFR = reshape(meanTFR,numFrequencies, windowLength);
        
        meanTFR = meanTFR(1:freqStop,timeStart:timeStop);
        
        
        
        % Calculate what theta windows contain slow and fast gamma episodes
        % Indices to theta cycles with gamma.
        [slowGammaWindowsInd, fastGammaWindowsInd, slowAndFastGammaWindowsInd] = getThetaWindowsWithGamma(thetaWindowStart, thetaWindowStop, slowGammaPeakIndNonOverlap, fastGammaPeakIndNonOverlap);
        
        onlySlowGammaWindowsInd = setdiff(slowGammaWindowsInd, slowAndFastGammaWindowsInd);
        slowGammaTFR = TFR(onlySlowGammaWindowsInd,:,:);
        
        meanSlowGammaTFR = mean(slowGammaTFR,1);
        meanSlowGammaTFR = reshape(meanSlowGammaTFR,numFrequencies, windowLength);
        meanSlowGammaTFR = meanSlowGammaTFR(1:freqStop,timeStart:timeStop);
        
        
        onlyFastGammaWindowsInd = setdiff(fastGammaWindowsInd, slowAndFastGammaWindowsInd);
        fastGammaTFR = TFR(onlyFastGammaWindowsInd,:,:);
        
        meanFastGammaTFR = mean(fastGammaTFR,1);
        meanFastGammaTFR = reshape(meanFastGammaTFR,numFrequencies, windowLength);
        meanFastGammaTFR = meanFastGammaTFR(1:freqStop,timeStart:timeStop);
        
        
        fastAndSlowGammaCombined = union(slowGammaWindowsInd, fastGammaWindowsInd);
        combinedGammaTFR = TFR(fastAndSlowGammaCombined,:,:);
        
        meanCombinedGammaTFR = mean(combinedGammaTFR,1);
        meanCombinedGammaTFR = reshape(meanCombinedGammaTFR,numFrequencies, windowLength);
        meanCombinedGammaTFR = meanCombinedGammaTFR(1:freqStop,timeStart:timeStop);
        
        disp('Calculating power of slow and fast gamma')
        % Calculate the mean power of the gamma
        nfft = 56;
        [meanPowerSlowGamma,freq] = doPwelch(mean(slowGammaWindowsEEG), eegData1{2}, nfft);
        [meanPowerFastGamma] = doPwelch(mean(fastGammaWindowsEEG), eegData1{2}, nfft);
        % Calculate the sum of the power
        powerSlowGamma = sum(meanPowerSlowGamma);
        powerFastGamma = sum(meanPowerFastGamma);
        % Normalize the power
        meanPowerSlowGamma = meanPowerSlowGamma / powerSlowGamma;
        meanPowerFastGamma = meanPowerFastGamma / powerFastGamma;
        [~,slowGammaStartInd] = min((freq - p.slowGammaFrequency(1)).^2);
        slowGammaStartInd = slowGammaStartInd(1);
        [~,slowGammaStopInd] = min((freq - p.slowGammaFrequency(2)).^2);
        slowGammaStopInd = slowGammaStopInd(1);
        [~,fastGammaStartInd] = min((freq - p.fastGammaFrequency(1)).^2);
        fastGammaStartInd = fastGammaStartInd(1);
        [~,fastGammaStopInd] = min((freq - p.fastGammaFrequency(2)).^2);
        fastGammaStopInd = fastGammaStopInd(1);
        % Axona recording system
        meanPowerSlowGammaEEG1SlowBand = meanPowerSlowGamma(slowGammaStartInd:slowGammaStopInd,:);
        meanPowerSlowGamma = mean(meanPowerSlowGammaEEG1SlowBand);
        meanPowerFastGammaEEG1FastBand = meanPowerFastGamma(fastGammaStartInd:fastGammaStopInd,:);
        meanPowerFastGamma = mean(meanPowerFastGammaEEG1FastBand);
        try
        Gamma_Results.Mean_power_of_slow_gamma(z,1) = meanPowerSlowGamma;
        Gamma_Results.Mean_power_of_fast_gamma(z,1) = meanPowerFastGamma;
        catch
            keyboard
        end
        disp('Calculating number of consecutive theta cycles with slow and fast gamma')
    

        [~, thetaWindowStart, thetaWindowStop, thetaIndices] = cut2thetaShort(eegData1{1},eegData1{2},0);
        [slowGammaWindowsInd, fastGammaWindowsInd] = getThetaWindowsWithGamma(thetaWindowStart, thetaWindowStop, slowGammaPeakIndNonOverlap,...
            fastGammaPeakIndNonOverlap);
        
        [thetaWithSlowGammaInRowHist, thetaWithFastGammaInRowHist] = thetaCyclesInRowWithGamma(slowGammaWindowsInd,...
            fastGammaWindowsInd, thetaIndices);
        
        % Total number of theta cycles in the recording
        numThetaCycles = length(thetaIndices);
        % Probability for a theta cycle to contain slow gamma
        pSlow = length(slowGammaWindowsInd) / numThetaCycles; 
        expectedThetaWithSlowGammaInRowHist = zeros(10,1);
        for ii = 1:10
            expectedThetaWithSlowGammaInRowHist(ii) = pSlow^ii * (1-pSlow);
        end
        
        % Probability for a theta cycle to contain fast gamma
        pFast = length(fastGammaWindowsInd) / numThetaCycles; 
        expectedThetaWithFastGammaInRowHist = zeros(10,1);
        for ii = 1:10
            expectedThetaWithFastGammaInRowHist(ii) = pFast^ii * (1-pFast);
        end
        try
        Gamma_Results.Probability_of_theta_containing_fast_gamma(z,1) = pFast;
        Gamma_Results.Probability_of_theta_containing_slow_gamma(z,1) = pSlow;
        
        Gamma_Results.Obs_minus_expected_probability_fast_gamma_on_theta_cycles_in_row{z} = thetaWithFastGammaInRowHist(1,:)' - expectedThetaWithFastGammaInRowHist(:,1);
        Gamma_Results.Obs_minus_expected_probability_slow_gamma_on_theta_cycels_in_row{z} = thetaWithSlowGammaInRowHist(1,:)' - expectedThetaWithSlowGammaInRowHist(:,1);
        catch
            keyboard
        end
        disp('Calculating correlation between slow and fast gamma')
     
        [~, thetaWindowStart, thetaWindowStop] = cut2thetaShort(eegData1{1},eegData1{2},0);

        % Calculate the slow fast gamma correlation and p-value
        [slowFastCorrcoeff,slowFastP] = slowFastGammaCorrelation(slowGammaPeakIndNonOverlap, fastGammaPeakIndNonOverlap, thetaWindowStart, thetaWindowStop);
        
        slowFastCorrcoeff = slowFastCorrcoeff(2);
        slowFastP = slowFastP(2);
        try
        Gamma_Results.Correlation_r_between_slow_and_fast_gamma(z,1) = slowFastCorrcoeff;
        Gamma_Results.Correlation_p_between_slow_and_fast_gamma(z,1) = slowFastP;
        catch
            keyboard
        end
       
        %% PLOTZ
        if plotfig == 1  
          % Plot the theta phase gamma episode histogram
        thetaPhaseBinSlowGammaDouble = [thetaPhaseBinSlowGamma; thetaPhaseBinSlowGamma];
        thetaPhaseBinFastGammaDouble = [thetaPhaseBinFastGamma; thetaPhaseBinFastGamma];
        plotAxis = 15:30:705;
     
        figure(1)
        bar(plotAxis,thetaPhaseBinSlowGammaDouble,'hist')
        title('Number of slow gamma episodes on different phases of theta (Normalized)')
        xlim([0, 720]);
        fName = sprintf('%s%s%s',File,'_slow_gamma_theta_phase_count');
        imageStore(figure(1),p.imageFormat,fName,300);
        
        figure(2)
        bar(plotAxis,thetaPhaseBinFastGammaDouble,'hist')
        title('Number of fast gamma episodes on different phases of theta (Normalized)')
        xlim([0, 720]);
        fName = sprintf('%s%s%s',File,'_fast_gamma_theta_phase_count');
        imageStore(figure(2),p.imageFormat,fName,300);
        
         % Make the figure of the CRF    
        figure(3)
        imagesc(freqVec1(1:300,:),freqVec2(:),CRF(:,1:300))
        axis xy
        colorbar
        xlabel('Frequency_p_h_a_s_e [Hz]')
        ylabel('Frequency_p_o_w_e_r [Hz]')
        colormap jet
        fName = sprintf('%s%s%s',File,'_CRF');
        imageStore(figure(3),p.imageFormat,fName,300);
        
         % This figure can be used to find the exact border values for the slow
        % and fast gamma. The intersection points between the mean line and the
        % graph would define the borders, but in som case the CRF is noisy and
        % the borders can not be found. In that case the predefine borders must
        % be used (slow gamma 25-55 Hz and fast gamma 65-140 Hz). Per default
        % the 25-55 Hz and 65-140 Hz bands are used.
        figure(4)
        clf
        plot(freqVec2, thetaCRF)
        hold on
        line([0, max(freqVec2)], [meanCoherence, meanCoherence],'color','k');
        hold off
        xlabel('Frequency [Hz]')
        ylabel('Coherence')
        fName = sprintf('%s%s%s',File,'_Thetaband_CRF');
        imageStore(figure(4),p.imageFormat,fName,300);
        Gamma_Results.freqvec{z} = freqVec2; 
        Gamma_Results.theta_crf{z} = thetaCRF;
        
        figure(5)
        bar(thetaWithSlowGammaInRowHist)
        title('Slow gamma');
        fName = sprintf('%s%s%s',File,'_thetaCyclesInRowWithSlowGamma');
        imageStore(figure(5),p.imageFormat,fName,300);
        
        Gamma_Results.theta_fast_in_row{z} = thetaWithFastGammaInRowHist;
        Gamma_Results.theta_slow_in_row{z} = thetaWithSlowGammaInRowHist;
        
        figure(6)
        bar(thetaWithFastGammaInRowHist)
        title('Fast gamma')
        fName = sprintf('%s%s%s',File,'_thetaCyclesInRowWithFastGamma');
        imageStore(figure(6),p.imageFormat,fName,300);
        close all
        else
        end
        process = strcat('Finished processing session','_',num2str(z),'_of_',num2str(height(EGF)));
        disp(process)
        Gamma_Results.File = EGF.path;
        clearvars -except  Gamma_Results z EGF p plotfig
        end
      
        save Gamma_Results_hires Gamma_Results