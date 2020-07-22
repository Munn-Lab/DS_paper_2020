% Calculates place map for cell spike firing during gamma episodes
function [map, xAxis, yAxis] = rateMapForFiringDuringGamma(posx, posy, post, cellTS, gammaPeakInd, spikeEegInd, Fs, p)

% Find size of arena
minX = nanmin(posx);
maxX = nanmax(posx);
minY = nanmin(posy);
maxY = nanmax(posy);
xLength = maxX - minX;
yLength = maxY - minY;

numSpikes = length(cellTS);

if numSpikes == 0
    % No spike for this cell, make an empty map and return
    numBins = round(xLength / p.binWidth);
    xAxis = zeros(numBins,1);
    for ii = 1:numBins
        xAxis(ii) = minX + (ii-1) * p.binWidth;
    end
    yAxis = xAxis;
    map = zeros(numBins);
    return
end


numGammaWindows = length(gammaPeakInd);
gammaSpikeTs = zeros(numSpikes,1);
offset = round(0.2 * Fs);
counter = 0;

% Locate spikes within the gamma windows
for ii = 1:numGammaWindows
    start = gammaPeakInd(ii) - offset;
    stop = gammaPeakInd(ii) + offset;
    
    ind = find(spikeEegInd >= start & spikeEegInd <= stop);
    if ~isempty(ind)
        N = length(ind);
        
        gammaSpikeTs(counter+1:counter+N) = cellTS(ind);
        
        % Remove the spikes we added from the original array to avoid that
        % they get added again.
        cellTS(ind) = [];
        spikeEegInd(ind) = [];
        
        % Increment the spike counter
        counter = counter + N;
    end
end

% Shorten the array
gammaSpikeTs = gammaSpikeTs(1:counter);

% Calculate the postion for the spikes
[spkx,spky] = spikePos(gammaSpikeTs,posx,posy,post);

[map, ~, xAxis, yAxis] = rateMap(posx,posy,spkx,spky,p.binWidth,p.binWidth,minX,xLength,minY,yLength,p.sampleTime,p);