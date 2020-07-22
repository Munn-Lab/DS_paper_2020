% Calculates instanteous theta filtered frequency
function EEG_Properties = eeg_amplitude_and_frequency(flt_eeg, peakInd, Fs)

numPeaks = length(peakInd);

% 1: Peak Index
% 2: Peak Time
% 3: Instantinous frequency
% 4: Instantinous amplitude
EEG_Properties = zeros(numPeaks,4);
EEG_Properties(:,1) = peakInd;
EEG_Properties(:,2) = peakInd / Fs;

for ii = 2:numPeaks
    % Instantaneous frequency
    EEG_Properties(ii,3) = 1 / (peakInd(ii)/Fs - peakInd(ii-1)/Fs);
    
    % Instantaneous amplitude
    EEG_Properties(ii,4) = nanmax(flt_eeg(peakInd(ii-1):peakInd(ii))) - nanmin(flt_eeg(peakInd(ii-1):peakInd(ii)));
end

