% Find windows of the EEG where there exist theta osciliation, short thetea
% window = 200 ms
% Laura Lee Colgin.
function [W, kstart_save, kstop_save, W_index] = cut2thetaShort(D,Fs,peak)

%ts_eeg = ts_eeg / 1000000;
Ftheta = 8;    % Hz
%Ftol = 0.005;  % s
Ftol = 0.025;  % s Modified June, 2009 to have less strict theta detection criterion, as requested by Nature reviewers
assert(0 <= peak);
assert(360 >= peak);
% Set the offset time to 100 ms
offsetTime = 0.1;

kstart = 1 - floor(offsetTime*Fs);         % I changed this to be only 1 theta cycle
kstop =  1 + floor(offsetTime*Fs);         %
time = (1:kstop-kstart+1)/Fs - offsetTime; % 

W = zeros(length(time),10000);  %for raw data

thetaBP = fftbandpass_gamma(D,Fs,Ftheta-5,Ftheta-4,Ftheta+4,Ftheta+5);    % bandpass for theta


W_index = [];
kstart = 0; 
j = 0;


offset = floor(offsetTime*Fs);

for k = 2:length(D)-1     % iterate through all samples of the EEG

    if (...
        ((225 <= peak && 315 >  peak) && (thetaBP(k-1) > 0 && thetaBP(k) < 0))... %downslope
    ||  (( 45 <= peak && 135 >  peak) && (thetaBP(k-1) < 0 && thetaBP(k) > 0))... %upslope
    ||  (( 45 >  peak || 315 <= peak) && ((thetaBP(k-1) > thetaBP(k))  && (thetaBP(k+1) > thetaBP(k))))... %trough  
    ||  ((135 <= peak && 225 >  peak) && ((thetaBP(k-1) < thetaBP(k))  && (thetaBP(k+1) < thetaBP(k))))... %peak
    )
        kstartold = kstart;
        % kstart is sample index of now - 100 ms
        kstart = k - offset;
        % kstop is sample index of now + 100 ms
        kstop =  k + offset;
        
        if kstart > 0 && kstop < length(thetaBP) - 512     % window is inside EEG
            % dk is number of samples since last theta cycle
            dk = kstart - kstartold;
            
            if dk < Fs/Ftheta + Fs*Ftol && dk > Fs/Ftheta - Fs*Ftol % if number of samples is within the region around (Fs*Ftol) what is expected for theta (Fs/FTheta)
                j = j + 1;                                          % then it's probably a theta cycle

                W(:,j) = D(kstart:kstop);       % and we make a new cut
                W_index(j) = k;
                
                kstart_save(j) = kstart;
                kstop_save(j) = kstop;
                
            end
        end
    end
end
W = W(:,1:j);
