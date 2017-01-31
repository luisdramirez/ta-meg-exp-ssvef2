%%% Creating Contrast Checkerboard Target Stimulus

%% add paths
addpath(genpath('/Users/luisramirez/Documents/CarrascoLabMEG/vistadisp'))
addpath('/Users/luisramirez/Documents/CarrascoLabMEG/ta-meg-exp/TAPilot')

%% 
displayName = 'Carrasco_L1'; % 'meg_lcd','Carrasco_L2','Carrasco_L1'
d = loadDisplayParams(displayName);

% pixelsPerDegree = 1/d.pixelSize;
pixelsPerCm = 1/d.pixelSize;
degPerCm = atan(1/d.distance)*180/pi;
pixelsPerDegree = pixelsPerCm/degPerCm;

screenWidth = d.numPixels(1); % (px)
screenHeight = d.numPixels(2); % (px)
cx = round(screenWidth/2);
cy = round(screenHeight/2);

%%
contrasts = [0.9 0.1 0.5]; 
phases = [0 pi];
pixelsPerDegree = 100;
stimSize = 3;
spatialFreq = 2;
orientation = 0;
blurRadius = 0.2;
backgroundColor = 128/255;

%% Make the stimuli
for iPhase = 1:numel(phases)
    for iContrast = 1:numel(contrasts)
        phase = phases(iPhase);
        contrast = contrasts(iContrast);
        

        c1 = buildColorGrating(pixelsPerDegree, [stimSize stimSize], ...
            spatialFreq, orientation, phase, contrast, ...
            1, 'bw', 1, 1);
        c2 = buildColorGrating(pixelsPerDegree, [stimSize stimSize], ...
            spatialFreq, orientation+90, phase, contrast, ...
            1, 'bw', 1, 1);
        if iPhase==1
            c = c1==c2;
        elseif iPhase==2
            c = c1~=c2;
        end
        s{iPhase, iContrast} = contrast*(double(c)-.5) + .5; %CONTRAST CHANGE HERE
        
        stim{iPhase, iContrast} = maskWithAnnulus(s{iPhase,iContrast}, ...
            length(s{iPhase,iContrast}), ...
            0, blurRadius, backgroundColor); %STIMULUS IMAGE GENERATED HERE
    end
end

% generate all combinations of phase, contrast, and side of screen
for iP1 = 1:numel(phases)
    for iC1 = 1:numel(contrasts)
        for iP2 = 1:numel(phases)
            for iC2 = 1:numel(contrasts)
                
                stimMatrix{iP2,iC2} = ... %iP1,iC1,
                    [stim{iP2,iC2}]; %stim{iP1,iC1} 
                
                stimIDs(iP2,iC2) = ... %iP1,iC1,
                    iP2*10 + iC2; % [phase left, contrast left, phase right, contrast right] iP1*1000 + iC1*100 +
            end
        end
    end
end

% put into images, keeping track of the IDs
for iIm = 1:numel(stimMatrix)
    images(:,:,iIm) = stimMatrix{iIm};
    imageIDs(iIm,1) = stimIDs(iIm);
end
imageIDHeaders = {'right-phase', 'right-contrast'}; %'left-phase', 'left-contrast', 

% add blank stimulus
images(:,:,end+1) = ones(size(images(:,:,1)))*backgroundColor;
imageIDs(end+1) = 0;

% stim n,n+1 are same contrast different phase
% stim n,n+2 are different contrast same phase
% each row is different phase, each column is different orientation
%% Create target stimulus

nPos = 8; %number of x,y pairs to place the gaussian
nCBPhases = 2; %number of checker board phases
nConds = 2; %number of contrast conditions 
nIms = nPos * nCBPhases * nConds; %total number of images to create 

xmax = size(im,1); ymax = size(im,2);

% Generate guassian center coordinates 
cx2 = 0; cy2 = 0;  %origin of coordiantes
r = xmax/4; %radius of center coordinates
theta = 22.5:45:360; theta = deg2rad(theta); 
x0 = cx2 + r * cos(theta); %generate x coordinates
y0 = cx2 + r * sin(theta); %generate y coordinates
gaussCoords = [x0' y0']; %store coordinates
shuffledCoords = gaussCoords(randperm(size(gaussCoords,1)),:); %chuffle coordinates

% 2D Gaussian parameters
gaussWidth = xmax;
gaussHeight = ymax;
sigma = 25; %radius of gaussian
gaussAmp = 1; %gaussian amplitude


% Generate empty cell arrays to store target images
circBlurs = cell(nPos,1);
maskedIms = cell(nPos, nConds, nCBPhases);

% target will be on for 4 frames
% make target images for 2 phases, increment/decrement, 8 locations

% Generate circular mask | use make2DGaussianCentered(w, h, x0, y0, sigma, gaussAmp)
for i=1:nPos
    circBlur = make2DGaussianCentered(gaussWidth, gaussHeight,  shuffledCoords(i,1), shuffledCoords(i,2), sigma, gaussAmp); 
    circBlurs{i} = circBlur;
%     figure(i)
%     imshow(circBlurs{i})
end

% Create target image for each phase and condition
for i= 1:nConds %for each condition (column in stim, page in maskedIms)
    for j= 1:nCBPhases %for each phase (row in stim, column in maskedIms)
        im = stim{j,i};
        for k= 1:nPos 
            maskedIm = cat(3, im, circBlurs{k}); %make masked image
            maskedIms{k,j,i} = maskedIm;
        end
    end
end
%% Display 
AssertOpenGL;
Screen('Preference', 'SkipSyncTests', 1);
screenNumber=max(Screen('Screens'));
x_win = 600;
y_win = x_win-200;
[w, wRect]=Screen('OpenWindow',screenNumber, 0,[0 0 cx cy],32,2);
Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% FIXATION
% make frame rect
rectDiam = 10;
baseRect = [0 0 rectDiam rectDiam];
centeredRect = CenterRectOnPointd(baseRect, cx/2, cy/2);
fixationColor = [0 0 0]*255;
fixationMaxDiameter = max(baseRect) * 1.01;

% make fixation rect
frameRect = [0 0 rectDiam*2 rectDiam*2];
centeredRect2 = CenterRectOnPointd(frameRect, cx/2, cy/2);
frameColor = [1 1 0]*255; %change outer ring color [1 1 0]=yellow, [1 0 0] = red
frameMaxDiameter = max(frameRect) *1.01;

% background stim/image
backgroundIms = cell(nCBPhases,1);
backgroundIms{1} = stim{5}; %'backgroundIm' is the fixed contrast. last column
backgroundIms{2} = stim{6};
for i= 1:nConds
    for j= 1:nCBPhases
        for k=1:nPos
            % Make target texture
            tex = Screen('MakeTexture', w, maskedIms{k,j,i}*255);
            % Make background texture
            bgtex = Screen('MakeTexture', w, backgroundIms{j}*255);

            % Blank sceen
            gray=GrayIndex(screenNumber); 
            Screen('FillRect',w, uint8(gray));
            Screen('Flip', w);

            % Show image
            Screen('DrawTexture', w, bgtex);
            Screen('DrawTexture', w, tex); % sRect, dRect);
            Screen('FrameOval', w, frameColor, centeredRect2, frameMaxDiameter)
            Screen('FillOval', w, fixationColor, centeredRect, fixationMaxDiameter)
            Screen('FrameOval', w, fixationColor, centeredRect2)
            Screen('Flip',w);
            pause(1)
        end
    end
end
sca