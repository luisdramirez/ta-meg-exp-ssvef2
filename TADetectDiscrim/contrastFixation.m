function fixation = contrastFixation(display, frame)
% Generates new dot fixation; called in 'showScanStimulus'

%7 = green; 6 = red;

% FIXATION
% make frame rect
rectDiam = 7;
baseRect = [0 0 rectDiam rectDiam];
centeredRect = CenterRectOnPointd(baseRect, display.fixX, display.fixY);


fixationColor = [0 0 0]*255;
fixationMaxDiameter = max(baseRect) * 1.01;

% make fixation rect
frameRect = [0 0 rectDiam*2 rectDiam*2];
centeredRect2 = CenterRectOnPointd(frameRect, display.fixX, display.fixY);

%change outer ring color
if frame == 5;
    frameColor = [1 1 1]*255;
elseif frame == 1
    frameColor = [0.5 0.5 0.5*255];
elseif frame == 6
    frameColor =  [1 0 0]*255;
elseif frame == 7
    frameColor = [0 1 0]*255;
elseif frame == 8
    frameColor = [0 0 1]*255;
end

frameMaxDiameter = max(frameRect) *1.01;

% AssertOpenGL;
% Screen('Preference', 'SkipSyncTests', 1);
% screenNumber=max(Screen('Screens'));        
% [w, wRect]=Screen('OpenWindow',screenNumber, 0,[0 0 x_win y_win],32,2);
% Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Blank sceen
% gray=GrayIndex(screenNumber); 
% Screen('FillRect',w, uint8(gray));

% w = 10;
Screen('FrameOval', display.windowPtr, frameColor, centeredRect2, frameMaxDiameter)
Screen('FillOval', display.windowPtr, fixationColor, centeredRect, fixationMaxDiameter)
Screen('FrameOval', display.windowPtr, fixationColor, centeredRect2)
return
% Screen('Flip', w);