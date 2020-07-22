function [information,sparsity,selectivity] = shuffle_func(Axona_Output,p) 
for j = 1:size(Axona_Output,2)
        if isempty(Axona_Output(j).cellTS)
            continue
        end
        cellTS = Axona_Output(j).cellTS; post = Axona_Output(j).post; posx = Axona_Output(j).posx;
        posx2 = Axona_Output(j).posx2; posy = Axona_Output(j).posy; posy2 = Axona_Output(j).posy2;
        eeg = Axona_Output(j).eeg; samplerate = Axona_Output(j).eeg_samplerate;
        spiketrain = computeFR(cellTS,post);
        File = strcat(Axona_Output(j).Filename,'_CELL_',num2str(Axona_Output(j).cellnum));
        speed = speed2D(posx,posy,post);
        bad_ind = speed < 2 | speed > 100;
        posx = posx(~bad_ind); posy = posy(~bad_ind);
        posx2 = posx2(~bad_ind); posy2 = posy2(~bad_ind); post_b = post(~bad_ind);       
        [spkx,spky,~] = spikePos(cellTS,posx,posy,post_b);      
        for ii = 1:length(spkx)
         a = min(spkx);
        b = max(spkx);
        spkx_shuff(ii,1) = (b-a).*rand + a;
        posx_shuff(ii,1) = (b-a).*rand + a;
        a = min(spky);
        b = max(spky);
        spky_shuff(ii,1) = (b-a).*rand + a;
        posy_shuff(ii,1) = (b-a).*rand + a;
        end
        xStart = min(posx_shuff); yStart = min(posy_shuff); xLength = max(posx_shuff) - min(posx_shuff); yLength = max(posy_shuff) - min(posy_shuff);
        sampleTime = post(2) - post(1);
        [map, posPdf, ~, ~] = ratemapAdaptiveSmoothing(posx_shuff, posy_shuff, spkx_shuff, spky_shuff, xStart, xLength, yStart, yLength, sampleTime, p, 1);
        [information(j,1),sparsity(j,1),selectivity(j,1)] = mapstat(map,posPdf);
end
end