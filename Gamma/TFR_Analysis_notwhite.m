function [thetaWindows, TFR, thetaWindowStart, thetaWindowStop] = TFR_Analysis_notwhite(EEG, freqVec, Fs, width)


[thetaWindows, thetaWindowStart, thetaWindowStop] = cut2theta(EEG,Fs,0);



% % Do the whitening of the signal for illustration purposes
% whiteW = zeros(size(thetaWindows));
% numWindows = size(thetaWindows,2);
% 
% for ii = 1:numWindows
%     whiteW(:,ii) = prewhitening(thetaWindows(:,ii),5);
% end
% 

% Calculate the TFR using Morlet wavelet
% TFR is 3-dimensional
% 1 dim = Number of theta windows
% 2 dim = Numger of frequencies to test
% 3 dim = Length of window in samples
TFR = traces2TFR(thetaWindows,freqVec,Fs,width);




%
%    Copyright (C) 2005 by Ole Jensen, F. C. Donders Centre for
%    for Cognitive Neuroimaging, PO Box 6500, Nijmegen, The 
%    Netherlands. Do not distribute without permission. 
%


%We used order = 5
% do this after you cut the theta cycles
% for picture of averaged theta cycles (i.e., Figure 1Ee, use the original,
% unwhitened theta cycles