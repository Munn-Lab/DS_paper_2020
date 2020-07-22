
% Calculates mean theta amplitude for speed bins
function [speedAxis,thetaFreqBin,thetaAmpBin] =  speedThetaBinning(speed,thetaFreq,thetaAmp, p)

% Maximum speed
maxSpeed = nanmax(speed);

numSpeedBins = maxSpeed / p.binWidthSpeed;

speedAxis = zeros(round(numSpeedBins),1);
thetaAmpBin = zeros(round(numSpeedBins),1);
thetaFreqBin = zeros(round(numSpeedBins),1);
start = 0;
stop = p.binWidthSpeed;

for ii = 1:numSpeedBins
    speedAxis(ii) = (start+stop) / 2;
    thetaAmpBin(ii) = nanmean(thetaAmp(speed >= start & speed < stop));
    thetaFreqBin(ii) = nanmean(thetaFreq(speed >= start & speed < stop));
    start = start + p.binWidthSpeed;
    stop = stop + p.binWidthSpeed;
end
end


