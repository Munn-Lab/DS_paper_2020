% Calculates a 2 dimensional rate map. The map is smoothed with a Gaussian
% smoothing kernel implemented with a boxcar lowpass filter, that effectively
% approximates a Gaussian filter.
%
% posx          x-coordinate for all the position samples in the recording
% posy          y-coordinate for all the position samples in the recording
% spkx          x-coordinate for all the spikes for a specific cell in the recording
% spky          y-coordinate for all the spikes for a specific cell in the recording
% xBinWidth     Bin width for the bins in map in the x-direction [cm]
% yBinWidth     Bin width for the bins in map in the Y-direction [cm] (Usually the same as the x bin width)
% xLength       Length of the arena in the x-direction [cm](for cylinder this equals the diameter)
% yLength       Length of the arena in the y-direction [cm] (for cylinder this equals the diameter)
% sampleTime    Sample duarion. For Axona it is 0.02 sec, for NeuraLynx it is 0.04 sec
%
% Version 1.0   
% 13. Dec 2007
%
% Version 1.1   Optimization for speed
% 20. Jan. 2012
%
% (c) Raymond Skjerpeng, Centre for the Biology of Memory, NTNU, 2007.
function [map, rawMap, xAxis, yAxis, timeMap] = rateMap(posx,posy,spkx,spky,xBinWidth,yBinWidth,xStart,xLength,yStart,yLength,sampleTime,p)
p.smoothingMode = 0;
% Number of bins in each direction of the map
numBinsX = ceil(xLength/xBinWidth);
numBinsY = ceil(yLength/yBinWidth);

% Allocate memory for the maps
spikeMap = zeros(numBinsY,numBinsX);
timeMap = zeros(numBinsY,numBinsX);

xAxis = zeros(numBinsX,1);
yAxis = zeros(numBinsY,1);



% Overall objective:
% Foreach (x-bin,y-bin) pair,
% count in spikemap/timemap (xbins x ybins) the number of places where
% (spky is in ybin and spkx is in xbin)
% (posy is in ybin and posx is in xbin)

% Bucketsort spikes and samples into regular bins
% Fortranesqe base1-indexing, add 1 for good measure
spkx_bin_idx = floor(((spkx - xStart) / xBinWidth)) + 1;
spky_bin_idx = floor(((spky - yStart) / yBinWidth)) + 1;
timex_bin_idx = floor(((posx - xStart) / xBinWidth)) + 1;
timey_bin_idx = floor(((posy - yStart) / yBinWidth)) + 1;
for n=1:length(spkx_bin_idx)
    ii = spkx_bin_idx(n);
    jj = spky_bin_idx(n);
    if ( ii>0 && ii<=numBinsX && jj>0 && jj<=numBinsY)
        spikeMap((numBinsY-jj+1),ii) = spikeMap((numBinsY-jj+1),ii) + 1;
    end
end
for n=1:length(timex_bin_idx)
    ii = timex_bin_idx(n);
    jj = timey_bin_idx(n);
    if ( ii>0 && ii<=numBinsX && jj>0 && jj<=numBinsY)
        timeMap((numBinsY-jj+1),ii) = timeMap((numBinsY-jj+1),ii) + 1;
    end
end

% Transform the number of spikes to time
timeMap = timeMap * sampleTime;

rawMap = spikeMap ./ timeMap;
rawMap(timeMap < p.minBinTime) = NaN;

if p.smoothingMode == 0
    % Smooth the spike and time map
    spikeMap = boxcarSmoothing(spikeMap);
    timeMap = boxcarSmoothing(timeMap);
else
    % Smooth the spike and time map
    spikeMap = boxcarSmoothing3x3(spikeMap);
    timeMap = boxcarSmoothing3x3(timeMap);
end


% Calculate the smoothed rate map
map = spikeMap ./ timeMap;

map(timeMap<p.minBinTime) = NaN;

% Set the axis
start = xStart + xBinWidth/2;
for ii = 1:numBinsX
    xAxis(ii) = start + (ii-1) * xBinWidth;
end
start = yStart + yBinWidth/2;
for ii = 1:numBinsY
    yAxis(ii) = start + (ii-1) * yBinWidth;
end



% Gaussian smoothing using a boxcar method
function sMap = boxcarSmoothing(map)

% Load the box template
box = boxcarTemplate2D();

% Using pos and phase naming for the bins originate from the first use of
% this function.
[numPhaseBins,numPosBins] = size(map);

sMap = zeros(numPhaseBins,numPosBins);

for ii = 1:numPhaseBins
    for jj = 1:numPosBins
        for k = 1:5
            % Phase index shift
            sii = k-3;
            % Phase index
            phaseInd = ii+sii;
            % Boundary check
            if phaseInd<1
                phaseInd = 1;
            end
            if phaseInd>numPhaseBins
                phaseInd = numPhaseBins;
            end
            
            for l = 1:5
                % Position index shift
                sjj = l-3;
                % Position index
                posInd = jj+sjj;
                % Boundary check
                if posInd<1
                    posInd = 1;
                end
                if posInd>numPosBins
                    posInd = numPosBins;
                end
                % Add to the smoothed rate for this bin
                sMap(ii,jj) = sMap(ii,jj) + map(phaseInd,posInd) * box(k,l);
            end
        end
    end
end


% Gaussian boxcar template
function box = boxcarTemplate2D()

% Gaussian boxcar template
box = [0.0025 0.0125 0.0200 0.0125 0.0025;...
       0.0125 0.0625 0.1000 0.0625 0.0125;...
       0.0200 0.1000 0.1600 0.1000 0.0200;...
       0.0125 0.0625 0.1000 0.0625 0.0125;...
       0.0025 0.0125 0.0200 0.0125 0.0025;];
   
   