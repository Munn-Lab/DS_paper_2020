
CELLS = table;
[file,path] = uigetfile('*.mat','Choose the combined output files','multiselect','on');
cd(path)
temp = strsplit(path,'/');
Animal = temp{length(temp)-1};
p_struct
for h = 1:length(file)
    CELLS.path(h,1) = strcat(string(path),file{h});
end
Results = struct;
All_Results = struct;
for z = 1:height(CELLS)
    load(CELLS.path(z,1));
    for j = 1:size(Axona_Output,2)
      
        cellTS = Axona_Output(j).cellTS;
          if isempty(cellTS)
            j = j + 1;
            cellTS = Axona_Output(j).cellTS;
        end
        
        post = Axona_Output(j).post; posx = Axona_Output(j).posx;
        posx2 = Axona_Output(j).posx2; posy = Axona_Output(j).posy; posy2 = Axona_Output(j).posy2;
        eeg = Axona_Output(j).eeg; samplerate = Axona_Output(j).eeg_samplerate;
        spiketrain = computeFR(cellTS,post);
        File = strcat(Axona_Output(j).Filename,'_CELL_',num2str(Axona_Output(j).cellnum));
        speed = speed2D(posx,posy,post);
        bad_ind = speed < 2 | speed > 100;
        posx = posx(~bad_ind); posy = posy(~bad_ind);
        posx2 = posx2(~bad_ind); posy2 = posy2(~bad_ind); post_b = post(~bad_ind);
        [spkx,spky,spkInd] = spikePos(cellTS,posx,posy,post_b);
        xStart = min(posx); yStart = min(posy); xLength = max(posx) - min(posx); yLength = max(posy) - min(posy);
        sampleTime = post(2) - post(1);
        [map, posPdf, rowAxis, colAxis] = ratemapAdaptiveSmoothing(posx, posy, spkx, spky, xStart, xLength, yStart, yLength, sampleTime, p, 1);
        [information,sparsity,selectivity] = mapstat(map,posPdf);
        theta = fftbandpass(eeg,samplerate,4,5,12,13);
        [~, peakInd] = thetaPhase(theta);
        thetaProperties = thetaFrequency(theta, peakInd, samplerate);
        [posTheta,posAmplitude]  = thetaFrequencyPosSample(post,thetaProperties);
        p.binWidthSpeed = 2;
        [speedAxis, thetaFreqBin, thetaAmpBin] =  speedThetaBinning(speed, posTheta, posAmplitude, p);
        speedAxis = speedAxis(1:50); thetaFreqBin = thetaFreqBin(1:50); thetaAmpBin = thetaAmpBin(1:50);
        [theta_output] = polyfit(speedAxis(~isnan(thetaFreqBin)),thetaFreqBin(~isnan(thetaFreqBin)),1); 
        theta_slope = theta_output(1);
        theta_intercept = theta_output(2);
        [deg_phase,~] = thetaPhase(theta); 
        phase = deg2rad(deg_phase); 
        singlet_pos = spiketrain == 1; doublet_pos = spiketrain == 2;  triplet_pos = spiketrain == 3;
        quadruplet_pos = spiketrain == 4;   more_pos = spiketrain > 4;
        num_singlets = sum(singlet_pos);    num_doublets = sum(doublet_pos);    num_triplets = sum(triplet_pos);
        num_quadruplets = sum(quadruplet_pos);  num_more_than_four = sum(more_pos);
        perc_single = num_singlets/length(cellTS);   perc_double = num_doublets/length(cellTS) * 100;
        perc_triple = num_triplets/length(cellTS) * 100;    perc_quad = num_quadruplets/length(cellTS) * 100;
        single_to_more_ratio = perc_single/(perc_double+perc_triple+perc_quad);
        phase_timestamp = (0:length(phase)-1)./250; 
        try
        for ii = 1:length(cellTS)
        diff = cellTS(ii) - phase_timestamp;
        [~,min_time(ii,1)] = min(abs(diff));
        end
        all_phases = phase(min_time);
        catch
            keyboard
        end
        phase_lock_vector = circ_r(all_phases(~isnan(all_phases))');
        if isnan(phase_lock_vector)
            keyboard
        end
        mean_phase = circ_mean((all_phases(~isnan(all_phases))'));
        mean_phase_degrees = rad2deg(mean_phase); 
   
        [Result_place] = place_cell_analysis_rm(posx,posy,post_b,cellTS,p);
  
        Results(j,1).Filename = strcat(Animal,'_',File);
        Results(j,1).mean_firing_rate = Result_place.mean_fr;
        Results(j,1).Number_of_fields = Result_place.number_of_fields;
        Results(j,1).mean_infield_rate = Result_place.mean_infield_rate_half_max;
        Results(j,1).mean_outfield_rate = Result_place.mean_outfield_rate_half_max;
        Results(j,1).place_field_percent_coverage = Result_place.percent_field_halfmax;
        Results(j,1).Information_content = information;
        Results(j,1).Sparsity = sparsity;
        Results(j,1).Selectivity = selectivity;
        Results(j,1).Theta_intercept = theta_intercept;
        Results(j,1).Theta_slope = theta_slope;
        Results(j,1).Phase_lock_vector = phase_lock_vector;
        Results(j,1).Mean_phase_direction = mean_phase_degrees;
        Results(j,1).Ratio_of_single_spikes = single_to_more_ratio;
        process = strcat('Finished processing cell','_',Results(j).Filename);
        disp(process)
%         name = strcat(Results(j).Filename,'_','_place_cell_results.mat');
%         save(name,Results);
        clear min_time all_phases phase_timestamp phase Result_place
    end
        clearvars -except Animal File All_Results CELLS p spiketrain z Results
%         filename = strcat(Results(j).Filename,'_Place_cell_Results.mat');
%         save(filename,'Results')
    if z == 1
        All_Results = Results;
    else
    All_Results = vertcat(All_Results,Results);
    end
    clear Results
end
