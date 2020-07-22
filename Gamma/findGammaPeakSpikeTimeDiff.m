% Made by Laura Colgin.
function [diff_time_ms] = findGammaPeakSpikeTimeDiff(eeg_epoch_index2,start2,stop2,non_overlap_peak_ind_unique,Fs)

%eeg_epoch_index2 is the output of "CA1_gamma_modulation_of_EC_CA3_spikes"
%start2, stop2, and non_overlap_peak_ind_unique are the output variables of
%"code_for_slow_fast_gamma_indices_windows3"



N = length(eeg_epoch_index2);
spikes_in_window = zeros(20*N,1);
corresponding_gamma_peak = zeros(20*N,1);
sampCounter = 0;



for ii = 1:length(start2)
    ind = find(eeg_epoch_index2 >= start2(ii) & eeg_epoch_index2 <= stop2(ii));
    N = length(ind);
    spikes_in_window(sampCounter+1:sampCounter+N) = eeg_epoch_index2(ind);
    corresponding_gamma_peak(sampCounter+1:sampCounter+N) = non_overlap_peak_ind_unique(ii);
    sampCounter = sampCounter + N;
end

spikes_in_window = spikes_in_window(1:sampCounter);
corresponding_gamma_peak = corresponding_gamma_peak(1:sampCounter);

samples_diff = spikes_in_window - corresponding_gamma_peak;
diff_time = samples_diff./Fs;
diff_time_ms = round(diff_time.*1000);