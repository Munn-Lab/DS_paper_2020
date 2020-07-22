function [TFR] = TFR_frequency_band(W,Fs,width, f1, f2)

freq(1,:) = [f1 f2];

ll = freq(1,1):2:freq(1,2);
N = length(ll);
temp = zeros(N,size(W,1));

n = 0;
for l = ll
    n = n + 1;
    temp(n,:) = energyvecForTFR(l,detrend(W),Fs,width);
end
TFR(1,1:size(W,1)) = mean(temp, 1);