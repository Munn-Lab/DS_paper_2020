load('/Users/rob/Documents/OneDrive - National University of Ireland, Galway/DOWN_SYNDROME_PROJECT/DATA/Path_to_Cells_and_EGF.mat');
Results = struct;
All_Results = struct;
for z = 1:height(CELLS)
    load(CELLS.CELLS(z,1));
    for j = 1:size(Axona_Output,2)
        cellTS = Axona_Output(j).cellTS; post = Axona_Output(j).post; posx = Axona_Output(j).posx;
        posx2 = Axona_Output(j).posx2; posy = Axona_Output(j).posy; posy2 = Axona_Output(j).posy2;
        Ripple_logic = CELLS.Ripple_Logic{j,:};
        spiketrain = computeFR(cellTS,post);
        egf_timestamp = (0:length(Ripple_logic)-1)./4800;
        
        for g = 1:length(cellTS) %%% here we're finding the closest value to 0 when we subtract the time of the spike (cellTS) from the sample of EEG (phase_timestamp)
        diff = cellTS(g) - egf_timestamp;
        [~,min_time(g)] = min(abs(diff)); %% This finds the value closest to 0 and returns an index and the value. Here we're supressing the value (by using '~') and getting the index
        end
        
        during_ripple = Ripple_logic(min_time);
        num_spikes_during_ripple = sum(during_ripple);
        percent_during_ripple = (num_spikes_during_ripple/length(cellTS)) * 100;
        File = strcat(Axona_Output(j).Filename,'_CELL_',num2str(Axona_Output(j).cellnum));
        speed = speed2D(posx,posy,post);
        [Y,Ty] = resample(double(Ripple_logic),250,4800);
        
        for h = 1:length(Y)
            if Y(h) ~= 0
                logic_rip(h,1) = 1;
            else
                logic_rip(h,1) = 0;
            end
        end
        
        logic_rip = logical(logic_rip);
        ripple_speed = speed(logic_rip);
        Results(j,1).mean_speed_during_ripples = nanmean(ripple_speed);
        non_ripple_speed = speed(~logic_rip);
        Results(j,1).mean_speed_outside_ripples = nanmean(non_ripple_speed);
        p_val = ranksum(ripple_speed,non_ripple_speed);
        if p_val < 0.05
            Results(j,1).Speed_sig_diff = 1;
        else
            Results(j,1).Speed_sig_diff = 0;
        end
        an = strsplit(CELLS.CELLS(z),'/');
        ani = strsplit(an{end},'Combined_Output');
        animal = ani{1};
        Sesh = sprintf('%s%s%s',animal,'Cellnumber_',num2str(j));
%         bad_ind = speed < 2 | speed > 100;
%         posx = posx(~bad_ind); posy = posy(~bad_ind);
%         posx2 = posx2(~bad_ind); posy2 = posy2(~bad_ind); post_b = post(~bad_ind);
        [spkx,spky,spkInd] = spikePos(cellTS,posx,posy,post);
        try
        ripple_spkx = spkx(during_ripple);
        ripple_spky = spky(during_ripple);
        catch
            keyboard
        end
        Results(j,1).Filename = Sesh;
        Results(j,1).num_spikes_during_ripple = num_spikes_during_ripple;
        Results(j,1).total_spikes = length(cellTS);
        Results(j,1).percent_spikes_during_ripple = percent_during_ripple;
        Results(j,1).posx = posx; Results(j,1).posy = posy;
        Results(j,1).spkx = spkx; Results(j,1).spky = spky;
        Results(j,1).ripple_spkx = ripple_spkx; Results(j,1).ripple_spky = ripple_spky;
   
    end
        process = strcat('Finished processing cell','_',Results(j).Filename);
        disp(process)
        filename = strcat(Results(j).Filename,'_SWR_Results.mat');
        save(filename,'Results')
        clearvars -except Animal File CELLS p j z Axona_Output
end
       clearvars -except CELLS Animal File p z
