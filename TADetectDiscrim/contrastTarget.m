function [backgroundIms, maskedIms, target] = contrastTarget(display, target)

%%% Creating Contrast Target Stimulus

%% Make the stimuli
phases = target.phases; 
pixelsPerDegree = target.pixelsPerDegree;
stimSize = target.stimSize;
spatialFreq = target.spatialFreq;
orientation = target.orientation;
blurRadius = target.blurRadius;
backgroundColor = target.backgroundColor;
%radialCB = target.radialCB;

%%
for iPhase = 1:numel(phases)
    for iContrast = 1:numel(target.contrast)
        phase = phases(iPhase);
        contrast = target.contrast(iContrast);
        switch target.stimType
            case 'checkerboard'
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
                targ_s{iPhase, iContrast} = contrast*(double(c)-.5) + .5; %CONTRAST CHANGE HERE
            case 'bullseye'
                targ_s{iPhase, iContrast} = CreateSpiral(display, stimSize, spatialFreq, phase, contrast)./2 + .5; 
            case 'radialcb'
                targ_s{iPhase, iContrast} = makeRadialCheckerboard(pixelsPerDegree, stimSize, phase, contrast, ...
                    radialCB.thetaCycles, radialCB.E, radialCB.A, radialCB.b);
            otherwise
                error('stimType not recognized')
        end
        
        targ_stim{iPhase, iContrast} = maskWithAnnulus(targ_s{iPhase,iContrast}, ...
            length(targ_s{iPhase,iContrast}), ...
            0, blurRadius, backgroundColor); %STIMULUS IMAGE GENERATED HERE
        
    end
end

% generate all combinations of phase, contrast, and side of screen
for iP1 = 1:numel(phases)
    for iC1 = 1:numel(target.contrast)
        for iP2 = 1:numel(phases)
            for iC2 = 1:numel(target.contrast)
                
                targ_stimMatrix{iP2,iC2} = ... %iP1,iC1,
                    [targ_stim{iP2,iC2}]; %stim{iP1,iC1} 
                
                targ_stimIDs(iP2,iC2) = ... %iP1,iC1,
                    iP2*10 + iC2; % [phase left, contrast left, phase right, contrast right] iP1*1000 + iC1*100 +
            end
        end
    end
end

% put into images, keeping track of the IDs
for iIm = 1:numel(targ_stimMatrix)
    targ_images(:,:,iIm) = targ_stimMatrix{iIm};
    targ_imageIDs(iIm,1) = targ_stimIDs(iIm);
end
imageIDHeaders = {'right-phase', 'right-contrast'}; %'left-phase', 'left-contrast', 

% add blank stimulus
targ_images(:,:,end+1) = ones(size(targ_images(:,:,1)))*backgroundColor;
targ_imageIDs(end+1) = 0;

% stim n,n+1 are same contrast different phase
% stim n,n+2 are different contrast same phase
% each row is different phase, each column is different orientation
%% Create target stimulus

nPos = 8; %number of x,y pairs to place the gaussian
nPhases = 2; %number of checker board phases
nConds = 2; %number of contrast conditions 
nIms = nPos * nPhases * nConds; %total number of images to create 

coords = target.coords;
xmax = size(targ_stim{1},1); ymax = size(targ_stim{1},2);

% 2D Gaussian parameters
gaussWidth = xmax;
gaussHeight = ymax;
sigma = 25; %radius of gaussian
gaussAmp = 1; %gaussian amplitude


% Generate empty cell arrays to store target images
circBlurs = cell(nPos,1);
maskedIms = cell(nPos, nPhases , nConds);

% target will be on for 4 frames
% make target images for 2 phases, increment/decrement, 8 locations

% Generate circular mask | use make2DGaussianCentered(w, h, x0, y0, sigma, gaussAmp)
for frame=1:nPos
    circBlur = make2DGaussianCentered(gaussWidth, gaussHeight,  coords(frame,1), coords(frame,2), sigma, gaussAmp); 
    circBlurs{frame} = circBlur;
%     figure(i)
%     imshow(circBlurs{i})
end

% Create target image for each phase and condition
for frame= 1:nConds %for each condition (column in targStim, page in maskedIms)
    for j= 1:nPhases %for each phase (row in targStim, column in maskedIms)
        im = targ_stim{j,frame};
        for k= 1:nPos 
            maskedIm = cat(3, im, circBlurs{k}); %make masked image
            maskedIms{k,j,frame} = maskedIm;
        end
    end
end

% background stim/image
backgroundIms = cell(nPhases,1);
backgroundIms{1} = target.stim{1}; %'backgroundIm' is the fixed contrast. last column
backgroundIms{2} = target.stim{2};

%% Display 

% AssertOpenGL;
% Screen('Preference', 'SkipSyncTests', 1);
% screenNumber=max(Screen('Screens'));        
% [w, wRect]=Screen('OpenWindow',screenNumber, 0,[0 0 x_win y_win],32,2);
% Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% for i= 1:nConds
%     for j= 1:nPhases
%         for k=1:nPos
%             % Make target texture
%             tex = Screen('MakeTexture', w, maskedIms{k,j,i}*255);
%             % Make background texture
%             bgtex = Screen('MakeTexture', w, backgroundIms{j}*255);
% 
%             % Blank sceen
%             gray=GrayIndex(screenNumber); 
%             Screen('FillRect',w, uint8(gray));
%             Screen('Flip', w);
% 
%             % Show image
%             Screen('DrawTexture', w, bgtex);
%             Screen('DrawTexture', w, tex); % sRect, dRect);
% 
%             Screen('Flip',w);
%             pause(0.5)
%         end
%     end
% end
% sca


% % TESTING
% AssertOpenGL;
% Screen('Preference', 'SkipSyncTests', 1);
% screenNumber=max(Screen('Screens'));        
% [w, wRect]=Screen('OpenWindow',screenNumber, 0,[0 0 x_win y_win],32,2);
% Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% 
% % Blank sceen
% gray=GrayIndex(screenNumber); 
% Screen('FillRect',w, uint8(gray));
% 
% %Make target texture
% tex = Screen('MakeTexture', w, maskedIms{1,1,1}*255);
% %Make background texture
% bgtex = Screen('MakeTexture', w, backgroundIms{1}*255);
% 
% %Show image
% Screen('DrawTexture', w, bgtex);
% Screen('DrawTexture', w, tex);
% Screen('Flip',w);
