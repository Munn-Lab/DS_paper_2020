function [map, posPdf, rowAxis, colAxis] = ratemapAdaptiveSmoothing(posx, posy, spkx, spky, xStart, xLength, yStart, yLength, sampleTime, p, shape)
% Calculates an adaptive smoothed rate map as described in "Skaggs et al
% 1996 - Theta Phase Precession in Hippocampal Neuronal Population and the
% Compression of Temporal Sequences"
%
% Input arguments
%
% posx          x-coordinate for all the position samples in the recording
% posy          y-coordinate for all the position samples in the recording
%
% spkx          x-coordinate for all the spikes for a specific cell in the 
%               recording
% spky          y-coordinate for all the spikes for a specific cell in the
%               recording
%
% xStart        Minimum x-coordinate for the path
%
% yStart        Minimum y-coordinate for the path
%
% xLength       Length of the arena in the x-direction [cm](for cylinder 
%               this equals the diameter)
% yLength       Length of the arena in the y-direction [cm] (for cylinder
%               this equals the diameter)
% sampleTime    Sample duarion. For Axona it is 0.02 sec, for NeuraLynx it
%               is 0.04 sec
%
% p             Parameter list with p.binWidth and p.alpha value
%
% shape         Shape of the box. Square box = 1, Cylinder = 2.
%
% Output variables
%
% map           The adaptive smoothed map
%
% posPdf        The position probability density function
%
%
% Version 1.0
% 13. Jan. 2010
%
% (c) Raymond Skjerpeng, KI/CBM, NTNU, 2010.




% Number of bins in each direction of the map
numColBins = ceil(xLength/p.binWidth);
numRowBins = ceil(yLength/p.binWidth);

rowAxis = zeros(numRowBins,1);
for ii = 1:numRowBins
    rowAxis(numRowBins-ii+1) = yStart+p.binWidth/2+(ii-1)*p.binWidth;
end
colAxis = zeros(numColBins, 1);
for ii = 1:numColBins
    colAxis(ii) = xStart+p.binWidth/2+(ii-1)*p.binWidth;
end

maxBins = max([numColBins, numRowBins]);

map = zeros(numRowBins, numColBins);
posPdf = zeros(numRowBins, numColBins);


binPosX = (xStart+p.binWidth/2);

if shape(1) == 1
    for ii = 1:numColBins

        binPosY = (yStart + p.binWidth/2);
        
        for jj = 1:numRowBins

            radius = maxBins * p.binWidth;
            % Number of samples inside the circle
            n = insideCircle(binPosX, binPosY, radius, posx, posy);
            % Number of spikes inside the circle
            s = insideCircle(binPosX, binPosY, radius, spkx, spky);

            if maxBins > p.alphaValue/(n*sqrt(s))         
                
                n = 0;
                s = 0;
                for r = 1:maxBins
                    % Set the current radius of the circle
                    radius = r * p.binWidth;
                    % Number of samples inside the circle
                    n = insideCircle(binPosX, binPosY, radius, posx, posy);
                    % Number of spikes inside the circle
                    s = insideCircle(binPosX, binPosY, radius, spkx, spky);

                    if r >= p.alphaValue/(n*sqrt(s))         
                        break;
                    end

                end
            end
            
    
            % Set the rate for this bin
            map(numRowBins-jj+1,ii) = s/(n*sampleTime);
            posPdf(numRowBins-jj+1,ii) = n*sampleTime;
            binPosY = binPosY + p.binWidth;
        end 

        binPosX = binPosX + p.binWidth;
    end

else
    for ii = 1:numColBins

        binPosY = (yStart + p.binWidth/2);
        for jj = 1:numRowBins
            currentPosition = sqrt(binPosX^2 + binPosY^2);
            if currentPosition > shape(2)/2
                map(numRowBins-jj+1,ii) = NaN;
                posPdf(numRowBins-jj+1,ii) = NaN;
            else
                n = 0;
                s = 0;
                for r = 1:maxBins
                    % Set the current radius of the circle
                    radius = r * p.binWidth;
                    % Number of samples inside the circle
                    n = insideCircle(binPosX, binPosY, radius, posx, posy);
                    % Number of spikes inside the circle
                    s = insideCircle(binPosX, binPosY, radius, spkx, spky);

                    if r >= p.alphaValue/(n*sqrt(s))         
                        break;
                    end

                end
                % Set the rate for this bin
                map(numRowBins-jj+1,ii) = s/(n*sampleTime);
                posPdf(numRowBins-jj+1,ii) = n*sampleTime;
                
            end
            binPosY = binPosY + p.binWidth;
        end 

        binPosX = binPosX + p.binWidth;
    end
end

map(posPdf<0.100) = NaN;
posPdf = posPdf / nansum(nansum(posPdf));
end

function n = insideCircle(cx, cy, radius, pointX, pointY)

dist = sqrt((pointX-cx).^2 + (pointY-cy).^2);
n = length(dist(dist <= radius));
end