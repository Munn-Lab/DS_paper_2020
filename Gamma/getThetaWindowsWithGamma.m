function [slowGammaWindows, fastGammaWindows, slowAndFastGammaWindows] = getThetaWindowsWithGamma(thetaWindowStart, thetaWindowStop, slowGammaPeakIndNoOverlap, fastGammaPeakIndNoOverlap)

numSlow = length(slowGammaPeakIndNoOverlap);
numFast = length(fastGammaPeakIndNoOverlap);

N = length(thetaWindowStart);
slowGammaWindows = zeros(N,1);
fastGammaWindows = zeros(N,1);

for w = 1:N
    for n = 1:numSlow
        if slowGammaPeakIndNoOverlap(n) > thetaWindowStart(w) && slowGammaPeakIndNoOverlap(n) < thetaWindowStop(w)
            slowGammaWindows(w) = 1;
        end
    end
    
    for n = 1:numFast
        if fastGammaPeakIndNoOverlap(n) > thetaWindowStart(w) && fastGammaPeakIndNoOverlap(n) < thetaWindowStop(w)
            fastGammaWindows(w) = 1;
        end
    end
end

slowGammaWindows = find(slowGammaWindows == 1);
fastGammaWindows = find(fastGammaWindows == 1);

slowAndFastGammaWindows = intersect(slowGammaWindows, fastGammaWindows);
