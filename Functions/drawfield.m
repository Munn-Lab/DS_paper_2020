% drawfield(map,xAxis,yAxis,cmap,maxrate)
%
% This function will calculate an RGB image from the rate
% map. We do not just call image(map) and caxis([0 maxrate]),
% as it would plot unvisted parts with the same colour code
% as 0 Hz firing rate. Instead we give unvisited bins
% their own colour (e.g. gray or white).
%
% Version 1.0       Written by Sturla Molden
%
% Version 1.1       Added possibility to have different number of values
% 24. Sep. 2008     in the x and y dimension of the rate map.
%
% Version 1.2       Possible to have all zero rate
% 09.Jul.2009
%
% (c) Sturla Molden and Raymond Skjerpeng
function drawfield(map,rowAxis,colAxis,cmap,maxrate)
   
[numRow,numCol] = size(map);
plotmap = ones(numRow,numCol,3);
if maxrate == 0
    for ii = 1:numRow
        for jj = 1:numCol
            if isnan(map(ii,jj))
                plotmap(numRow-ii+1,jj,1) = 1; % give the unvisited bins a gray colour
                plotmap(numRow-ii+1,jj,2) = 1;
                plotmap(numRow-ii+1,jj,3) = 1;
            else
                 plotmap(numRow-ii+1,jj,1) = 0; % give the unvisited bins a gray colour
                 plotmap(numRow-ii+1,jj,2) = 0;
                 plotmap(numRow-ii+1,jj,3) = 0.5625;
            end
        end
    end
else
   for ii = 1:numRow
        for jj = 1:numCol
            if isnan(map(ii,jj))
                plotmap(numRow-ii+1,jj,1) = 1; % give the unvisited bins a gray colour
                plotmap(numRow-ii+1,jj,2) = 1;
                plotmap(numRow-ii+1,jj,3) = 1;
            else
                rgb = pixelcolour(map(ii,jj),maxrate,cmap);
                plotmap(numRow-ii+1,jj,1) = rgb(1);
                plotmap(numRow-ii+1,jj,2) = rgb(2);
                plotmap(numRow-ii+1,jj,3) = rgb(3);
            end
        end
    end
end
image(rowAxis,colAxis,plotmap);
set(gca,'YDir','Normal');
axis('image')


   
% This function calculates a colour for each bin
% in the rate map.
function rgb = pixelcolour(map,maxrate,cmap)

cmap1 = ...
[    0         0    0.5625; ...
     0         0    0.6875; ...
     0         0    0.8125; ...
     0         0    0.9375; ...
     0    0.0625    1.0000; ...
     0    0.1875    1.0000; ...
     0    0.3125    1.0000; ...
     0    0.4375    1.0000; ...
     0    0.5625    1.0000; ...
     0    0.6875    1.0000; ...
     0    0.8125    1.0000; ...
     0    0.9375    1.0000; ...
0.0625    1.0000    1.0000; ...
0.1875    1.0000    0.8750; ...
0.3125    1.0000    0.7500; ...
0.4375    1.0000    0.6250; ...
0.5625    1.0000    0.5000; ...
0.6875    1.0000    0.3750; ...
0.8125    1.0000    0.2500; ...
0.9375    1.0000    0.1250; ...
1.0000    1.0000         0; ...
1.0000    0.8750         0; ...
1.0000    0.7500         0; ...
1.0000    0.6250         0; ...
1.0000    0.5000         0; ...
1.0000    0.3750         0; ...
1.0000    0.2500         0; ...
1.0000    0.1250         0; ...
1.0000         0         0; ...
0.8750         0         0; ...
0.7500         0         0; ...
0.6250         0         0 ];

cmap2 = ...
[0.0417         0         0; ...
0.1250         0         0; ...
0.2083         0         0; ...
0.2917         0         0; ...
0.3750         0         0; ...
0.4583         0         0; ...
0.5417         0         0; ...
0.6250         0         0; ...
0.7083         0         0; ...
0.7917         0         0; ...
0.8750         0         0; ...
0.9583         0         0; ...
1.0000    0.0417         0; ...
1.0000    0.1250         0; ...
1.0000    0.2083         0; ...
1.0000    0.2917         0; ...
1.0000    0.3750         0; ...
1.0000    0.4583         0; ...
1.0000    0.5417         0; ...
1.0000    0.6250         0; ...
1.0000    0.7083         0; ...
1.0000    0.7917         0; ...
1.0000    0.8750         0; ...
1.0000    0.9583         0; ...
1.0000    1.0000    0.0625; ...
1.0000    1.0000    0.1875; ...
1.0000    1.0000    0.3125; ...
1.0000    1.0000    0.4375; ...
1.0000    1.0000    0.5625; ...
1.0000    1.0000    0.6875; ...
1.0000    1.0000    0.8125; ...
1.0000    1.0000    0.9375];

if strcmp(cmap,'jet')
  steps = (31*(map/maxrate))+1;
  steps = round(steps);
  if steps > 32
      steps = 32;
  end 
  rgb = cmap1(steps,:);
else
  steps = (31*(map/maxrate))+1;
  steps = round(steps);
  rgb = cmap2(steps,:);
end
