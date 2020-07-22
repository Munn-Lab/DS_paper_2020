function [CRF,freqVec1,freqVec2] = crossfreqCoh(S1,S2,freqVec,Fs,p)

width = p.waveletWidth;


% The length of the fft should be about 4*Fs
N = Fs * 2;
% For speeding up the fft we go to the next power of 2 for the length.
nfft = 2^nextpow2(N);
if nfft/N > 1.9
    nfft = 2^nextpow2(N-1);
end

freqVec2 = freqVec';
% Number of frequency to calculate the coherence for
N = length(freqVec);
% Calculate the length of the mscohere output
if mod(nfft,2) == 0
    N2 = nfft/2 + 1;
else
    N2 = (nfft+1)/2;
end
% Allocate memory for the CRF
CRF = zeros(N,N2);

E1 = energyvec(freqVec(1),S1',Fs,width);
[CRF(1,:),freqVec1] = mscohere(E1',S2,hanning(nfft),nfft/2,nfft,Fs);

parfor k=2:length(freqVec)
    E1 = energyvec(freqVec(k),S1',Fs,width);
    CRF(k,:)= mscohere(E1',S2,hanning(nfft),nfft/2,nfft,Fs);
end

CRF = CRF(:,1:floor(end/2));
freqVec1 = freqVec1(1:floor(end/2));

