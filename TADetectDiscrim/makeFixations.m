function [fixDot, fixFrame] = makeFixations(display)

% Generates new dot fixation; called in 'showScanStimulus'

%7 = green; 6 = red;
% pixelsPerDegree = 1/d.pixelSize;
pixelsPerCm = 1/display.pixelSize;
degPerCm = atan(1/display.distance)*180/pi;
pixelsPerDegree = pixelsPerCm/degPerCm;

% FIXATION
% make frame rect
rectDiam = stimulus.fixDiam*pixelsPerDegree; %4 = 0.15;
baseRect = [0 0 rectDiam rectDiam];
centeredRect = CenterRectOnPointd(baseRect, display.fixX, display.fixY);


% Fixation and Frame colors; fixSeq (1,5,6,7,8,9)
frameColor = [0 0 0]*255; % black 
fixationColors = 255*[1 1 1 % white (5)
    0.5 0.5 0.5 % gray (1, blank)
    1 0 0 % red (6, incorrect) 
    0 1 0 % green (7, correct)
    0 0 1 % blue (8, no response)
    0.3 0.3 0.3]; % (9)

fixationMaxDiameter = max(baseRect)*1.01;

% make fixation rect
frameRect = [0 0 rectDiam*2 rectDiam*2];
centeredRect2 = CenterRectOnPointd(frameRect, display.fixX, display.fixY);
frameMaxDiameter = max(frameRect)*1.01;

%Generate fixation dots and store

fixDot = Screen('FrameOval', display.windowPtr, fixationColor, centeredRect2, frameMaxDiameter);
% Generation fixation frame and store
fixFrame = Screen('FrameOval', display.windowPtr, frameColor, centeredRect2);
return