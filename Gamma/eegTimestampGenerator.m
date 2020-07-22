function newEegData = eegTimestampGenerator(eegData)


    % EEG is from Axona system
        system = 0;
        newEegData{1} = eegData{1,1};

% Number of EEG samples
numSamples = length(eegData{1});
% Sampling rate
Fs = 250;
    newEegData{2} = Fs;
    % The Axona system has no timestamps for the EEG signal, but the
    % sampling rate is constant and can be used to make timestamps.
    maxTime = (numSamples-1) / Fs;
    eegTS = 0:1/Fs:maxTime;
    newEegData{4} = eegTS;
    % Axona data will always be converted to volts in the databaseMaker
    % program.
    newEegData{5} = 1;
