function [peak_ind_unique, non_overlap_peak_ind_unique, start2, stop2, non_overlap_gamma_windows_EEG1] = gammaEpisodes_one_channel(EEG, f1, f2, Fs)

%EEG = ch3_rev_begin1;

%f1 = 100;    %  frequency band of interest
%f2 = 140;

numEegSamples1 = length(EEG);

[TFR] = TFR_frequency_band(EEG,Fs,5, f1, f2);


TFR_z = zscore(TFR);
high = find(TFR_z > 2);
bp = fftbandpass_gamma(EEG,Fs,f1-2,f1,f2,f2+2);

%   length of the window of interest (time * sampling frequency)
detect_length = round(0.16 * Fs);
% start of the window of interest
detect_start = round(0.08 * Fs);

N = size(high,2);
peak_ind = zeros(N,1);
offset = round(0.2 * Fs);

kc = 0;
for k = 1:N %   for all samples with a high power
    % Determine start of window of interest
    start = high(k) - detect_start;
    start2 = high(k) - offset;
    % Determine stop of window of interest
    stop = start + detect_length;               
    stop2 = high(k) + offset;
    if start2 > 0 && stop2 <= numEegSamples1  % window has to be inside the whole EEG
        % Counter for windows of interest inside the whole EEG
        kc = kc + 1;
        % Determine the index of the maximal value in the bandpassed EEG
        [~, max_ind]= max(bp(start:stop)); 
        % Add the offset, in peak_ind are now all the sample indices with high power
        peak_ind(kc) = start + max_ind - 1;
    end
end
peak_ind = peak_ind(1:kc);

% We don't want duplicates
peak_ind_unique = unique(peak_ind);

N = length(peak_ind_unique);



non_overlap_peak_ind_unique = zeros(N,1);
counter = 0;
for n=2:length(peak_ind_unique)
    if peak_ind_unique(n-1) + (round(Fs/10)) < peak_ind_unique(n); %we require also that the peaks be separated by 100 ms
        counter = counter + 1;
        non_overlap_peak_ind_unique(counter) = peak_ind_unique(n-1);
    end
end

non_overlap_peak_ind_unique = non_overlap_peak_ind_unique(1:counter);
offset = round(0.2 * Fs);
numWindowSamples = 2 * offset + 1;
non_overlap_gamma_windows_EEG1 = zeros(counter,numWindowSamples);
non_overlap_gamma_windows_EEG2 = zeros(counter,numWindowSamples);
start2 = zeros(counter,1);
stop2 = zeros(counter,1);

for k = 1:counter

    % Start and stop sample for this window
    start = non_overlap_peak_ind_unique(k) - offset;
    stop = non_overlap_peak_ind_unique(k) + offset;
    
    if start > 0 && stop > 0
        

        % EEG samples for the window
        non_overlap_gamma_windows_EEG1(k,1:numWindowSamples) = EEG(start:stop);
        non_overlap_gamma_windows_EEG2(k,1:numWindowSamples) = EEG2(start:stop);

        start2(k) = start;
        stop2(k) = stop;
    end
end
