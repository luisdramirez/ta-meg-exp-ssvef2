function contrastFixation(display, stimulus, frame)
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


frameColor = [0 0 0]*255;
fixationMaxDiameter = max(baseRect)*1.01;

% make fixation rect
frameRect = [0 0 rectDiam*2 rectDiam*2];
centeredRect2 = CenterRectOnPointd(frameRect, display.fixX, display.fixY);

%change outer ring color
if stimulus.fixSeq(frame) == 5
    fixationColor = [1 1 1]*255;
elseif stimulus.fixSeq(frame) == 1
    fixationColor = [0.5 0.5 0.5]*255;
elseif stimulus.fixSeq(frame) == 6
    fixationColor =  [1 0 0]*255;
elseif stimulus.fixSeq(frame) == 7
    fixationColor = [0 1 0]*255;
elseif stimulus.fixSeq(frame) == 8
    fixationColor = [0 0 1]*255;
elseif stimulus.fixSeq(frame) == 9
    fixationColor = [0.3 0.3 0.3]*255;
end

frameMaxDiameter = max(frameRect)*1.01;

% AssertOpenGL;
% Screen('Preference', 'SkipSyncTests', 1);
% screenNumber=max(Screen('Screens'));        
% [w, wRect]=Screen('OpenWindow',screenNumber, 0,[0 0 x_win y_win],32,2);
% Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Blank sceen
% gray=GrayIndex(screenNumber); 
% Screen('FillRect',w, uint8(gray));

% w = 10;
Screen('FrameOval', display.windowPtr, fixationColor, centeredRect2, frameMaxDiameter)
% Screen('FillOval', display.windowPtr, fixationColor, centeredRect, fixationMaxDiameter)
Screen('FrameOval', display.windowPtr, frameColor, centeredRect2)
return
% Screen('Flip', w);