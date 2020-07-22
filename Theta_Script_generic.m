function [Results] = Theta_Script_generic

CELLS = table;
[file,path] = uigetfile('*.mat','Choose the combined output files','multiselect','on');
cd(path)
temp = strsplit(path,'/');
Animal = temp{length(temp)-1};
p_struct
    if ~iscell(file)
        file_t{1} = file;
        clear file
        file = file_t;
        clear file_t
    end
for h = 1:length(file)
    CELLS.path(h,1) = strcat(string(path),file{h});
end

All_Results = struct;
for z = 1:height(CELLS)
load(CELLS.path(z));
for j = 1:size(Axona_Output,2)
Results = struct;
%% Calculate theta phase of EEG

lfp = Axona_Output.eeg;
samplerate = Axona_Output.eeg_samplerate;
cellTS = Axona_Output(j).cellTS;
theta = fftbandpass(lfp,samplerate,4,5,12,13);
[deg_phase,peak_positions] = thetaPhase(theta);
% thetaProperties = thetaFrequency(theta, peak_positions,samplerate);
phase = deg2rad(deg_phase);
phase_timestamp = (0:length(phase)-1)./samplerate;


for g = 1:length(cellTS) %%% here we're finding the closest value to 0 when we subtract the time of the spike (cellTS) from the sample of EEG (phase_timestamp)
diff = cellTS(g) - phase_timestamp;
[~,min_time(g)] = min(abs(diff)); %% This finds the value closest to 0 and returns an index and the value. Here we're supressing the value (by using '~') and getting the index
end

spike_phases = phase(min_time);
spike_amplitude = theta(min_time);

name{1,:} = strcat(Axona_Output(j).Filename,'_','CELL_',num2str(Axona_Output(j).cellnum),'_of_',num2str(size(Axona_Output,2)));
mean_phase = circ_mean(spike_phases');
spike_lfp_phase_lock_vector = circ_r(spike_phases');
[circ_p,circ_z] = circ_rtest(spike_phases');
mean_amplitude = nanmean(abs(spike_amplitude));

Results.Session(j,1) = name;
Results.Number_of_spikes(j,1) = length(cellTS);

Results.theta_phase(j,1) = rad2deg(mean_phase+pi);
Results.phase_lock(j,1) = spike_lfp_phase_lock_vector;
Results.amplitude_at_spikes(j,1) = mean_amplitude;
Results.rao_Z(j,1) = circ_z; Results.rao_pval(j,1) = circ_p;
keyboard
% clear cellTS lfp samplerate theta min_time diff
end

di = strcat('Done_session_',num2str(z),'_of_',num2str(height(CELLS)));
disp(di);
if z == 1
    All_Results = Results;
else
    try
All_Results = vertcat(All_Results,Results);
    catch
        keyboard
    end
end

clear Results

end

