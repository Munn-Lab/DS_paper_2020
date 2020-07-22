% Parameters for axona gui
fprintf('%s','P Struct Loading');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% head_direction_and_gridness.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Size in centimeters for the bins in the ratemap
p.binWidth = 2.5;

% Method for smoothing the 2-D rate maps
% 0 = Gaussian smoothing with boxcar filter
% 1 = Adaptive smoothing
p.smoothingMethodRateMaps = 1;
p.hdSmoothingMode = 1;

% Size of the smoothing window when using flat smoothing window for the
% head direction map. The size is the total span of the filter. Please make
% the size an odd integer multippel of the head direction bin width
% (p.hdBinWidth). If this is not the case the program will round the number
% of to make it a odd multippel of the head direction bin width,
p.hdSmoothingWindowSize = 14.5; % [degrees]
p.hdBinWidth = 3; % [degrees]
p.hdBinWidth2 = 3; % [degrees]

% Threshold value for the orientation calculation (Made by Tor and Hanne)
p.gridOrientationThreshold = 0.5;

% A cell must have a gridness score at least as high as this for it to
% qualify.
p.gridnessScoreThreshold = 0.34;

% The mean vector length for a cell must be at least this high for the
% cell to qualify as a head direction cell.
p.meanVectorLengthThreshold = 0.20;

% Sets if the angle beween gridness lines should be used to pass cells as
% gridness cells or not.
% 1 = Use angle span
% 0 = Do not use it. Cells pass only based on gridness score.
p.gridnessAngleCriteriaMode = 0;

% Angle span when seeing if the grid orientation is good enough
p.angleSpanThreshold = 15;

% Minimum radius used in the auto-correlogram when finding the best
% gridness score
p.minRadius = 20; % [cm]

% Increment step of the radius when calculating what radius gives the best
% gridness score. 
p.radiusStep = p.binWidth; % [cm]

% When calculating gridness score based on the best radius, the program
% calculates the gridness score as an average over adjacent radii. This
% parameter sets the number of radii to use. The number must be odd. The
% middle radius will be set to be the best radius.
p.numGridnessRadii = 3;

% Threshold value used when defining the centre field. Border of centre
% field is defined as the place where the correlation falls under this
% threshold or the correlation start to increase again.
p.correlationThresholdForCentreField = 0.2;

% Mode = 0: Gridness is calculated as the mean correlation at 60 and 120
%           degrees minus the mean correlation at 30, 90 and 150 degrees.
% Mode = 1: Gridness is calculated as the minimum correlation at 60 and 120
%           degrees minus the maximum correlation at 30, 90 and 150
%           degrees.
p.gridnessCalculationMode = 1;

% Minimum allowed width of the correlogram disk. I.e the distance from the
% centre radius to the maximum radius
p.minDiskWidth = 20; % [cm]

% Alpha value for the adaptive smoothing rate map calculation. In use for
% the spatial information calculation
p.alphaValue = 10000;

% Same as p.alphaValue, but for the head direction map
p.hdAlphaValue = 10000;

% 0 =   Use direction of path as the head direction. Direction is only
%       calculated from 1 diode.
% 1 =   Use both LED to calculate the head direction.
p.headDirectionMode = 1;
p.conjunctify = 0;
p.minBinTime = 0.020; % [sec]

% 0 = Mean firing direction
% 1 = Peak firing direction
p.preferedDirecetionMode = 0;
p.hdMapLineWidth = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theta_analyser_rob.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bin width for the autocorrelogram
p.autocorrelogramBinWidth = 0.002; % [Sec]
% Length of the autocorrelogram. 
p.diagramLength = 0.500; % [Sec]
% Span of the power spectrum figure.
p.powerSpectrumSpan = 20; % [Hz]
% Minimum number of spikes. If the number of spikes is lower the analysis
% is skipped
p.minimumNumSpikes = 100;
% Frequencies below this value will not be used in calculating the mean
% power of the auto-correlogram.
p.lowestFrequency = 0; % [Hz]
% Frequency border values for the theta band
p.thetaLow = 4; % [Hz]
p.thetaHigh = 11; % [Hz]
% Low speed threshold. Segments of the path where the rat moves slower
% than this threshold will be removed. Set it to zero (0) to keep
% everything. Value in centimeters per second.
p.lowSpeedThreshold = 2.5; % [cm/s]
% High speed threshold. Segments of the path where the rat moves faster
% than this threshols will be removed. Set it to zero (0) to keep
% everything. Value in centimeters per second.
p.highSpeedThreshold = 100; % [cm/s]
% When speed filtering the eeg signal, signal discontinuities are
% introduced. This parameter sets the minimum length of a continous
% segment. Each discontinuity might introduce noise (spectral leakage) into
% the frequency spectre, so a low number of discontinuities are preferable.
p.minContinousSegmentLength = 0.5; % [sec]
% Picture storing format. 
% format = 1 -> bmp (24 bit)
% format = 2 -> png
% format = 3 -> eps
% format = 4 -> jpg
% format = 5 -> ill (Adobe Illustrator)
% format = 6 -> tiff (24 bit)
p.imageFormat = 2;
% Image folder
p.imageFolder = 'thetaImages';
% Bin width for the speed. 
p.binWidthSpeed = 5; % [cm/s]
% It is possible to only analyse part of the recording by setting these
% values. If both are set to zero the whole recording will be used. When
% one or both are set to non-zero values only data within the interval is
% used for analysing.
p.startTime = 0; % [second]              
p.stopTime = 0; % [second]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% place_cell_analysis.m               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.binWidth = 2.5;
% Method for smoothing the 2-D rate maps
% 0 = Gaussian smoothing with boxcar filter
% 1 = Adaptive smoothing
p.smoothingMethodRateMaps = 1;
p.hdSmoothingMode = 1;
% Size of the smoothing window when using flat smoothing window for the
% head direction map. The size is the total span of the filter. Please make
% the size an odd integer multippel of the head direction bin width
% (p.hdBinWidth). If this is not the case the program will round the number
% of to make it a odd multippel of the head direction bin width,
p.hdSmoothingWindowSize = 14.5; % [degrees]
p.hdBinWidth = 3; % [degrees]
p.hdBinWidth2 = 3; % [degrees]
% Threshold value for the orientation calculation (Made by Tor and Hanne)
p.gridOrientationThreshold = 0.5;
% A cell must have a gridness score at least as high as this for it to
% qualify.
p.gridnessScoreThreshold = 0.34;
% The mean vector length for a cell must be at least this high for the
% cell to qualify as a head direction cell.
p.meanVectorLengthThreshold = 0.20;
% Sets if the angle beween gridness lines should be used to pass cells as
% gridness cells or not.
% 1 = Use angle span
% 0 = Do not use it. Cells pass only based on gridness score.
p.gridnessAngleCriteriaMode = 0;
% Angle span when seeing if the grid orientation is good enough
p.angleSpanThreshold = 15;
% Minimum radius used in the auto-correlogram when finding the best
% gridness score
p.minRadius = 20; % [cm]
% Increment step of the radius when calculating what radius gives the best
% gridness score. 
p.radiusStep = p.binWidth; % [cm]
% When calculating gridness score based on the best radius, the program
% calculates the gridness score as an average over adjacent radii. This
% parameter sets the number of radii to use. The number must be odd. The
% middle radius will be set to be the best radius.
p.numGridnessRadii = 3;
% Threshold value used when defining the centre field. Border of centre
% field is defined as the place where the correlation falls under this
% threshold or the correlation start to increase again.
p.correlationThresholdForCentreField = 0.2;
% Mode = 0: Gridness is calculated as the mean correlation at 60 and 120
%           degrees minus the mean correlation at 30, 90 and 150 degrees.
% Mode = 1: Gridness is calculated as the minimum correlation at 60 and 120
%           degrees minus the maximum correlation at 30, 90 and 150
%           degrees.
p.gridnessCalculationMode = 1;
% Minimum allowed width of the correlogram disk. I.e the distance from the
% centre radius to the maximum radius
p.minDiskWidth = 20; % [cm]
% Alpha value for the adaptive smoothing rate map calculation. In use for
% the spatial information calculation
p.alphaValue = 10000;
% Same as p.alphaValue, but for the head direction map
p.hdAlphaValue = 10000;
% 0 =   Use direction of path as the head direction. Direction is only
%       calculated from 1 diode.
% 1 =   Use both LED to calculate the head direction.
p.headDirectionMode = 1;
p.conjunctify = 0;
p.minBinTime = 0.020; % [sec]
% 0 = Mean firing direction
% 1 = Peak firing direction
p.preferedDirecetionMode = 0;
p.hdMapLineWidth = 1;
p.lowestFieldRate = 0.5;
p.minNumBins = 2;
p.sampleTime = 0.02;
p.videoSamplingRate = 50;
p.percentile = 50; % [%]
% save('default_params.mat','p');