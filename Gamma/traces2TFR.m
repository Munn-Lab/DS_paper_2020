function [TFR,timeVec,freqVec] = traces2TFR(S,freqVec,Fs,width)
% function [TFR,timeVec,freqVec] = traces2TFR(S,freqVec,Fs,width);
%
% Calculates the average of a time-frequency energy representation of
% multiple trials using a Morlet wavelet method.
%
% Input
% -----
% S    : signals = time x Trials      
% freqVec    : frequencies over which to calculate TF energy        
% Fs   : sampling frequency
% width: number of cycles in wavelet (> 5 advisable)  
%
% Output
% ------
% t    : time
% f    : frequency
% B    : phase-locking factor = frequency x time
%     
%------------------------------------------------------------------------
% Ole Jensen, Brain Resarch Unit, Low Temperature Laboratory,
% Helsinki University of Technology, 02015 HUT, Finland,
% Report bugs to ojensen@neuro.hut.fi
%------------------------------------------------------------------------

%    Copyright (C) 2000 by Ole Jensen 
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You can find a copy of the GNU General Public License
%    along with this package (4DToolbox); if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


S = S';

numFrequencies = length(freqVec);
windowLength = size(S,2);
numWindows = size(S,1);

timeVec = (1:windowLength)/Fs;  


TFR = zeros(numWindows , numFrequencies, windowLength);
for ii = 1:numWindows
    B = zeros(numFrequencies,windowLength); 
    for jj = 1:length(freqVec)
        B(jj,:) = energyvecForTFR(freqVec(jj),detrend(S(ii,:)),Fs,width);
    end
    TFR(ii,:,:) = B;
end