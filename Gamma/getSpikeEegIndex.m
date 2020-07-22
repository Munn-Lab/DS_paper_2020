function spikeEegInd = getSpikeEegIndex(spikeTs, eegTs)

% TS_spikes is the time stamp for the spikes, output of the 'loadSpikes'
% function.
% ts_eeg is the time stamp for the EEG recordings (CSC file), one time
% stamp for every 512 points.  This variable is part of the output of the
% 'load_eeg' function.
% EEG is the CSC file containing the EEG recording to which the spikes will be aligned 


N = length(spikeTs);

if N == 0
    spikeEegInd = [];
    return
end

spikeEegInd = zeros(N,1);


for ii = 1:N
    [~, ind] = min((eegTs-spikeTs(ii)).^2);
    spikeEegInd(ii) = ind(1);
end
