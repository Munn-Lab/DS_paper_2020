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
Results = struct;
All_Results = struct;
for z = 1:height(CELLS)
    load(CELLS.path(z,1));
    for j = 1:size(Axona_Output,2)
        cellTS = Axona_Output(j).cellTS; post = Axona_Output(j).post; posx = Axona_Output(j).posx;
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
        [Result_place] = place_cell_analysis_rm(posx,posy,post_b,cellTS,p,File);
        Results(j,1).Filename = strcat(Animal,'_',File);
        Results(j,1).mean_firing_rate = Result_place.mean_fr;
        Results(j,1).Number_of_fields = Result_place.number_of_fields;
        Results(j,1).mean_infield_rate = Result_place.mean_infield_rate;
        Results(j,1).mean_outfield_rate = Result_place.mean_outfield_rate;
        Results(j,1).mean_infield_rate_half_max = Result_place.mean_infield_rate_half_max;
        Results(j,1).mean_outfield_rate_half_max = Result_place.mean_outfield_rate_half_max;
        Results(j,1).place_field_percent_half_max= Result_place.percent_field_halfmax;
   
    end
        process = strcat('Finished processing cell','_',Results(j).Filename);
        disp(process)
        filename = strcat(Results(j).Filename,'_Place_Cell_Results.mat');
        save(filename,'Results')
        clearvars -except Animal File CELLS p j z Axona_Output
end
       clearvars -except CELLS Animal File p z
