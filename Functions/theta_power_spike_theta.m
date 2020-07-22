[file,path] = uigetfile('.mat','multiselect','on');
cd(path);

if ~iscell(file)
    fi{1,1} = file;
    clear file
    file = fi;
    clear fi
end

for z = 73:size(file,2)
load(file{z});
for g = 1:size(Axona_Output,2)
    eeg = Axona_Output(g).hires_eeg;
    fs = Axona_Output(g).hires_samplerate;
    posx = Axona_Output(g).posx;
    posy = Axona_Output(g).posy;
    post = Axona_Output(g).post;
    cellTS = Axona_Output(g).cellTS;
    [spkx,spky,spkInd] = spikePos(cellTS,posx,posy,post);
    p_struct2;
%     [map, rawMap, xAxis, yAxis, timeMap] = rateMap(posx,posy,spkx,spky,2,2,min(posx),abs(min(posx))+max(posx),min(posy),abs(min(posy))+max(posy),0.02,p)
% 
%        [map, posPDF, rowAxis, colAxis]  = ratemapAdaptiveSmoothing(posx, posy, spkx, spky, min(posx),...
%            abs(min(posx))+max(posx),min(posy),abs(min(posy))+max(posy), 0.02, p, 1);
% stats = regionprops(map);
%     binary_map = imextendedmax(map,8);
%      [theta_eeg] = bandpass(eeg,6,12,fs);
     mean_power_in_theta_band = bandpower(eeg,fs,[6 12]);
     [thetaidx,p_val,mean_spike_freq,~,theta_index_ci,~,~,~,~] = theta_index(cellTS);
     wOutputs_theta = wFilter(eeg,fs,5,12,0);
     theta_eeg = wOutputs_theta{1};
     
%      wOutputs_logamma = wFilter(eeg,fs,25,55,0);
%      logamma_eeg = wOutputs_logamma{1};
%      
%      wOutputs_higamma = wFilter(eeg,fs,65,140,0);
%      higamma_eeg = wOutputs_higamma{1};
%      
%      eeg_timestamps = transpose(linspace(0,length(higamma_eeg)/fs,length(higamma_eeg)));
fi_name =  strcat(file{z},'_Cell_',num2str(g));
Results(g).File = fi_name;
Results(g).theta_power = mean_power_in_theta_band;
Results(g).spike_theta_index = thetaidx;
Results(g).spike_theta_pval = p_val;
Results(g).spike_theta_ci = theta_index_ci;
Results(g).mean_spike_freq = mean_spike_freq;

fil = strsplit(Results(g).File,'.');
filename = strcat(fil{1},'_cell_',num2str(g),'_Theta_spike_Results.mat');
save(filename,'Results');
name = strcat('Finished cell_',num2str(g),'_of_',num2str(size(Axona_Output,2)),'_of_session_',num2str(z),'_of_',num2str(size(file,2)));     
disp(name);
clearvars posx posy post cellTS eeg fs spkx spky spkInd
end
clearvars Results Axona_Output
end