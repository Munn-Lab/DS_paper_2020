% Calculates the phase bin histogram
function [phaseBin] = phaseBinHistogram(sPhase)

% Convert sPhase to degrees
sPhase = (sPhase+pi) * 360/(2*pi);

start = 0;
stop = 30;
phaseBin = zeros(12,1);
for ii = 1:12
    phaseBin(ii) = length(find(sPhase >= start & sPhase < stop));
    start = start + 30;
    stop = stop + 30;
end
