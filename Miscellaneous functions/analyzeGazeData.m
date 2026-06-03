function gazeStats = analyzeGazeData(version, elevation, varargin)
% ANALYZEGAZEDATA  Analyzes gaze data to detect and characterize saccades and fixations
%
%   gazeStats = analyzeGazeData(version, elevation)
%   Analyzes gaze position data to detect and characterize saccadic and fixation
%   events based on velocity thresholds.
%
%   Syntax
%   ------
%   gazeStats = analyzeGazeData(version, elevation)
%   gazeStats = analyzeGazeData(version, elevation, velocity)
%   gazeStats = analyzeGazeData(version, elevation, velocity, ...
%                               'binDuration', value, ...
%                               'thresholds', [upper lower])
%
%   Inputs
%   ------
%   version     - Vector of horizontal gaze positions (degrees)
%   elevation   - Vector of vertical gaze positions (degrees)
%   velocity    - (Optional) Vector of gaze velocities
%   binDuration - (Name-value) Bin duration = 1/sampling_rate (seconds)
%                 Default: 1
%   thresholds  - (Name-value) [upperThresh lowerThresh] in degrees/sec
%                 Default: [2000 25]
%
%   Outputs
%   -------
%   gazeStats   - Structure containing detected saccades and fixations with
%                 timing, position, and velocity information
%
%   Notes
%   -----
%   Saccades are detected when gaze velocity exceeds lowerThresh but stays
%   below upperThresh. All other samples are classified as fixations. Output
%   structure contains detailed metrics for each event including duration,
%   amplitude, mean velocity, and peak velocity.
%
%
%   Please cite
%   -----------
%     Hadjidimitrakis, K., Vaccari, F. E., De Vitis, M., Filippini, M., Diomedi, S., & Fattori, P. (2026).
%     Spontaneous oculomotor behavior sharpens eye position signals in parietal cortex.
%     [Manuscript under review].

    %% --- Check if velocity was provided ---
    if ~isempty(varargin) && isnumeric(varargin{1}) && ...
            length(varargin{1}) == length(version)
        
        velocity = varargin{1};
        varargin(1) = [];   % remove velocity from varargin
        velocityProvided = true;
    else
        velocityProvided = false;
    end

    %% --- Parse name-value arguments ---
    p = inputParser;
    addParameter(p, 'binDuration', 1, @isnumeric);
    addParameter(p, 'thresholds', [2000 25], @isnumeric);
    parse(p, varargin{:});

    binDuration = p.Results.binDuration;
    thresholds  = p.Results.thresholds;

    upperThresh = thresholds(1);
    lowerThresh = thresholds(2);

    %% --- Compute velocity only if not provided ---
    if ~velocityProvided
        velocity = sqrt(central_diff(version).^2 + ...
                        central_diff(elevation).^2);
    end

    velocityDegSec = velocity / binDuration;

    %% --- Initialize output structure ---
    gazeStats = struct();
    gazeStats.saccades = struct('startIdx', {}, 'endIdx', {}, 'duration', {}, ...
        'amplitudeX', {}, 'amplitudeY', {}, 'amplitudeXY', {}, ...
        'meanVelocity', {}, 'peakVelocity', {});
    gazeStats.fixations = struct('startIdx', {}, 'endIdx', {}, 'duration', {}, ...
        'amplitudeX', {}, 'amplitudeY', {}, 'amplitudeXY', {}, ...
        'meanVelocity', {}, 'peakVelocity', {});

    saccadeCount = 0;
    fixationCount = 0;
    i = 1;

    %% --- Main loop ---
    while i <= length(velocityDegSec)

        if velocityDegSec(i) >= lowerThresh && ...
           velocityDegSec(i) <= upperThresh

            % -------- Saccade --------
            startIdx = i;

            endIdx = startIdx;
            while endIdx < length(velocityDegSec) && ...
                  velocityDegSec(endIdx+1) >= lowerThresh && ...
                  velocityDegSec(endIdx+1) <= upperThresh
                endIdx = endIdx + 1;
            end

            saccadeCount = saccadeCount + 1;

            duration = (endIdx - startIdx) * binDuration;
            amplitudeX = version(endIdx) - version(startIdx);
            amplitudeY = elevation(endIdx) - elevation(startIdx);
            amplitudeXY = hypot(amplitudeX, amplitudeY);
            meanVelocity = mean(velocityDegSec(startIdx:endIdx));
            peakVelocity = max(velocityDegSec(startIdx:endIdx));

            gazeStats.saccades(saccadeCount) = struct( ...
                'startIdx', startIdx, ...
                'endIdx', endIdx, ...
                'duration', duration, ...
                'amplitudeX', amplitudeX, ...
                'amplitudeY', amplitudeY, ...
                'amplitudeXY', amplitudeXY, ...
                'meanVelocity', meanVelocity, ...
                'peakVelocity', peakVelocity);

            i = endIdx + 1;

        else

            % -------- Fixation --------
            startIdx = i;

            endIdx = startIdx;
            while endIdx < length(velocityDegSec) && ...
                 ~(velocityDegSec(endIdx+1) >= lowerThresh && ...
                   velocityDegSec(endIdx+1) <= upperThresh)
                endIdx = endIdx + 1;
            end

            fixationCount = fixationCount + 1;

            duration = (endIdx - startIdx + 1) * binDuration;
            amplitudeX = version(endIdx) - version(startIdx);
            amplitudeY = elevation(endIdx) - elevation(startIdx);
            amplitudeXY = hypot(amplitudeX, amplitudeY);
            meanVelocity = mean(velocityDegSec(startIdx:endIdx));
            peakVelocity = max(velocityDegSec(startIdx:endIdx));

            gazeStats.fixations(fixationCount) = struct( ...
                'startIdx', startIdx, ...
                'endIdx', endIdx, ...
                'duration', duration, ...
                'amplitudeX', amplitudeX, ...
                'amplitudeY', amplitudeY, ...
                'amplitudeXY', amplitudeXY, ...
                'meanVelocity', meanVelocity, ...
                'peakVelocity', peakVelocity);

            i = endIdx + 1;
        end
    end

    %% --- Metadata ---
    gazeStats.parameters.binDuration = binDuration;
    gazeStats.parameters.thresholds = thresholds;
    gazeStats.parameters.lowerThreshold = lowerThresh;
    gazeStats.parameters.upperThreshold = upperThresh;
    gazeStats.parameters.velocityProvided = velocityProvided;

    gazeStats.dataInfo.timeSeriesDuration = ...
        length(velocityDegSec) * binDuration;


end
