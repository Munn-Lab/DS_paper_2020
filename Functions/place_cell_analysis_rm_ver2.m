function [Results] = place_cell_analysis_rm(x,y,t,ts,p)

Results = struct;
fprintf('%s%s\n','Start analysing at ', datestr(now));

% Calculate the spike positions

[spkx,spky,spkInd] = spikePos(ts,x,y,t);

% Calculate the rate map
maxX = nanmax(x);
maxY = nanmax(y);
xStart = nanmin(x);
yStart = nanmin(y);
xLength = maxX - xStart + p.binWidth*2;
yLength = maxY - yStart +  p.binWidth*2;
    % Rate map with adaptive smoothing
    [map, posPDF, rowAxis, colAxis]  = ratemapAdaptiveSmoothing(x, y, spkx, spky, xStart, xLength, yStart, yLength, 0.02, p, 1);
    %[map, rawMap, rowAxis, colAxis, timeMap] = rateMap(x,y,spkx,spky,xBinWidth,yBinWidth,xStart,xLength,yStart,yLength,0.02,p)
    nn = length(rowAxis);
    
binary_map = imextendedmax(map,8);
meanrate = nansum(nansum( map .* posPDF ));
meansquarerate = nansum(nansum( (map.^2) .* posPDF ));
if meansquarerate == 0
   sparsity = NaN;
else
sparsity = meanrate^2 / meansquarerate;
end
maxrate = max(max(map));
if meanrate == 0
   selectivity = NaN;
else
   selectivity = maxrate/meanrate;
end
%% FIELD PROPERTIES

infield_bins = map(binary_map);
outfield_bins = map(~binary_map);

mean_infield_rate_imregionalmax = mean(infield_bins);
mean_outfield_rate_imregionalmax = mean(outfield_bins);
percent_field_imregionalmax = (length(infield_bins)/(size(map,1)*size(map,2))*100);
labelled_map = bwlabel(binary_map);
num_fields = max(max(labelled_map(:,:)));

for z = 1:num_fields
    bins_occ = find(labelled_map == z);
    field_size(z,1) = length(bins_occ)/(size(map,1)*size(map,2))*100;
end

%         avgRate = nanmean(nanmean(map(binsRow,binsCol)));
% %         avgRate_out = nanmean(nanmean(map(~binsRow,~binsCol)));
%         % Peak rate in field
%         peakRate = nanmax(nanmax(map(binsRow,binsCol)));
%         peakRate_out = nanmax(nanmax(setdiff(find(map),binsRow),setdiff(find(map),binsCol)));
        % Size of field

 keyboard
        Results.mean_fr = meanrate; Results.sparsity = sparsity;
        Results.number_of_fields = num_fields; Results.mean_infield_rate = mean_infield_rate_imregionalmax;
        Results.mean_outfield_rate = mean_outfield_rate_imregionalmax;
        Results.percent_field_coverage = percent_field_imregionalmax;
        Results.field_size_percent = field_size;
        Results.field_size = field_size; 

    fprintf('%s%s\n','Done analysing at ', datestr(now));
end
    

