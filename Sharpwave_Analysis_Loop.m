
EGF = table;
[file,path] = uigetfile('*.egf','Choose the combined output files','multiselect','on');
cd(path)
temp = strsplit(path,'/');
Animal = temp{length(temp)-1};
p_struct
Results = struct;
if isstring(file)
    for h = 1:size(file,1)
    EGF.path(h,1) = strcat(string(path),file);
    end
else
for h = 1:length(file)
    EGF.path(h,1) = strcat(string(path),file{h});
end
end

for z = 1:height(EGF)
    [status,eeg,fs,bytesPerSample] = readEGF(EGF.path(z,1));
    eegData{1,1} = eeg; eegData{1,2} = fs;
    plotfig = 0;
    [SWR, numSWR, freqSWR, fltrdEEG, ripple_locs] = detectSWR_KennethKay2016(eegData,plotfig);
    Results(z).ripples = SWR; Results(z).Number_ripples = numSWR;
    Results(z).ripples_per_sec = freqSWR; Results(z).ripple_logic =  ripple_locs;
    Results(z).cell = EGF.path(z,1);
    clearvars eegData eeg fs SWR ripple_locs
    process = strcat('Iteration number ','_',num2str(z),'_of_',num2str(height(EGF)));
    disp(process)
end

save Ripple_Results.mat Results

for t = 1:height(EGF)
        [status,eeg,fs,bytesPerSample] = readEGF(EGF.path(t,1));
        fltrdEEG = bandpass(eeg, 150,250, fs);
        eeg_ripples = fltrdEEG(Results(t).ripple_logic);
    [phase,peakPos] = eeg_phase(eeg_ripples);
    EEG_Properties = eeg_amplitude_and_frequency(eeg_ripples, peakPos, fs);
    Results(t).ripple_frequency = nanmean(EEG_Properties(:,3));
    Results(t).ripple_amplitude = nanmean(EEG_Properties(:,4));
%     TimeScale = 0:fs:length(eeg)
%     [x, y, z] = TimeFreq_Analysis(eeg, TimeScale, fs, 1)
end