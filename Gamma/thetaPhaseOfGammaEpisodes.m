% From Laura Colgin
function [phaseBin_theta_slow_gamma_norm, phaseBin_theta_fast_gamma_norm] = thetaPhaseOfGammaEpisodes(EEG, non_overlap_peak_ind_unique_slow, non_overlap_peak_ind_unique_fast, Fs, p)

% EEG_theta_filter = bandpass_gamma(EEG, p.thetaFrequency(1), p.thetaFrequency(2), Fs);
EEG_theta_filter = fftbandpass(EEG,Fs,p.thetaFrequency(1)-1,p.thetaFrequency(1),p.thetaFrequency(2)-1,p.thetaFrequency(2)); 
DTAS = hilbert(EEG_theta_filter);
theta_phase = angle(DTAS);


theta_phase_non_overlap_peak_ind_unique_slow = theta_phase(non_overlap_peak_ind_unique_slow);
theta_phase_non_overlap_peak_ind_unique_fast = theta_phase(non_overlap_peak_ind_unique_fast);




[phaseBin_theta_slow_gamma] = phaseBinHistogram(theta_phase_non_overlap_peak_ind_unique_slow);
phaseBin_theta_slow_gamma_norm = phaseBin_theta_slow_gamma/length(non_overlap_peak_ind_unique_slow);

[phaseBin_theta_fast_gamma] = phaseBinHistogram(theta_phase_non_overlap_peak_ind_unique_fast);
phaseBin_theta_fast_gamma_norm = phaseBin_theta_fast_gamma/length(non_overlap_peak_ind_unique_fast);
