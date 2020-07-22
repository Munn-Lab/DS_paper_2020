function [spkx,spky,spkInd] = spikePos(ts,posx,posy,post)
N = length(ts);
spkx = zeros(N,1);
spky = zeros(N,1);
spkInd = zeros(N,1);

count = 0;
for ii = 1:N
    tdiff = (post-ts(ii)).^2;
    [~,ind] = min(tdiff);

    % Check if spike is in legal time sone
    if ~isnan(posx(ind(1)))
        count = count + 1;
        spkx(count) = posx(ind(1));
        spky(count) = posy(ind(1));
        spkInd(count) = ind(1);
    end
end
spkx = spkx(1:count);
spky = spky(1:count);
spkInd = spkInd(1:count);
end