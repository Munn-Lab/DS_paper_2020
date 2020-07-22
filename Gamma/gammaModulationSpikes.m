function [eeg_epoch_index2_gamma_epochs_unique, eeg_epochs_grouped_ca1_gamma_epochs, sPhase2_ECCA3_cell_ca1_gamma_epochs, phaseBin_sPhase2_ECCA3_cell_ca1_gamma_epochs, phaseBin2_sPhase2_ECCA3_cell_ca1_gamma_epochs] = gammaModulationSpikes(eeg_CA1, peak_ind_unique, f1, f2, Fs, spikeEegInd )



N = length(peak_ind_unique);
eeg_epoch_index2_gamma_epochs = zeros(N,1);
counter2 = 0;
offset = round(0.2 * Fs);
for k = 1:N
    %   start of the windows that become averaged
    start = peak_ind_unique(k) - offset;
    %   stop of the windows that become averaged
    stop = start + 2 * offset;
    if 0 < start && length(eeg_CA1) >= stop %   windows have to be inside the whole EEG


        for m = 1:length(spikeEegInd)
            if spikeEegInd(m) > start && spikeEegInd(m) < stop 
                counter2 = counter2 + 1;
                eeg_epoch_index2_gamma_epochs(counter2) = spikeEegInd(m);
            end
        end
    end

end

eeg_epoch_index2_gamma_epochs = eeg_epoch_index2_gamma_epochs(1:counter2);

eeg_epoch_index2_gamma_epochs_unique = unique(eeg_epoch_index2_gamma_epochs);

[eeg_epochs_grouped_ca1_gamma_epochs, sPhase2_ECCA3_cell_ca1_gamma_epochs] = spike_triggered_windows_NLX_gamma40_phasevec_rob(eeg_epoch_index2_gamma_epochs_unique, eeg_CA1, f1, f2, Fs);


[phaseBin_sPhase2_ECCA3_cell_ca1_gamma_epochs, phaseBin2_sPhase2_ECCA3_cell_ca1_gamma_epochs] = make_spike_time_histogram_rep(sPhase2_ECCA3_cell_ca1_gamma_epochs);