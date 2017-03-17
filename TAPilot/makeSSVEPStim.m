function makeSSVEPStim(run)

%% run setup
% run = 3;
saveStim = 1;
saveFigs = 1;

%% file i/o
stimDir = 'stimuli';
stimFile = sprintf('taPilot%d', run);

%% screen setup
displayName = 'meg_lcd';
d = loadDisplayParams(displayName);
pixelsPerDegree = 1/d.pixelSize;
screenWidth = d.numPixels(1); % (px)
screenHeight = d.numPixels(2); % (px)
cx = round(screenWidth/2);
cy = round(screenHeight/2);

%% timing setup
refrate = 60; % (Hz)
blockDur = 7; % (s)
targetDur = 3/refrate; % (s)
targetLeadTime = 1; % (s) % no targets in first part of block
targetEndTime = 1; % (s) % no targets in last part of block
targetCushion = 1; % (s) % min interval between targets
maxTargetsPerBlock = 3;
attCueLeadTime = 0.5; % (s)
respDur = 1; % (s)
if refrate==75
    % 75 Hz SSVEP unit sequences: 4 frames (75/4=18.75 Hz) and 5 frames (75/5=15 Hz)
    fastUnit = [1 1 2 2]; % gives the phase (1 or 2) of each frame
    slowUnit = [1 1 1 2 2];
else
    % 60 Hz SSVEP unit sequences: 3 frames (60/3=20 Hz) and 4 frames (60/4=15 Hz)
    fastUnit = [1 2 2]; % gives the phase (1 or 2) of each frame
    slowUnit = [1 1 2 2];
end

%% blocks setup
blockNames = {'blank','fast-left','slow-left'};
% blockOrder = [1 2 1 2 1 3 1 3 1];
blockOrder = [1 2 1 2 1 3 1 3 1 2 1 2 1 3 1 3 1];
attBlockNames = {'no-att','att-left','att-right'};
% attBlockOrder = [1 2 1 3 1 2 1 3 1];
attBlockOrder = [1 2 1 3 1 2 1 3 1 2 1 3 1 2 1 3 1];
nBlocks = numel(blockOrder);

%% stim setup
stimSize = 8;
spatialFreq = 1;
orientation = 0;
stimContrast = 0.64;
targetContrast = 1;
contrasts = [stimContrast targetContrast];
blurRadius = 0.2;
backgroundColor = 128/255;
phases = [0 pi];

stimPos = [6 4]; % [x y]
stimSpacerWidth = (stimPos(1)-stimSize/2)*2;

%% trigger setup
triggerOption = 'conditionID'; % 'conditionID','combinatorial'

%% Store all stim params
% sorry this is kind of a horrible way to do this
p = v2struct(...
    displayName, pixelsPerDegree, screenWidth, screenHeight, cx, cy, ...
    refrate, blockDur, targetDur, targetLeadTime, targetEndTime, targetCushion, ...
    maxTargetsPerBlock, attCueLeadTime, respDur, fastUnit, slowUnit, ...
    blockNames, blockOrder, attBlockNames, attBlockOrder, ...
    stimSize, spatialFreq, orientation, stimContrast, targetContrast, ...
    contrasts, blurRadius, backgroundColor, phases, triggerOption);

%% Make the stimuli
for iPhase = 1:numel(phases)
    for iContrast = 1:numel(contrasts)
        phase = phases(iPhase);
        contrast = contrasts(iContrast);
        
        s{iPhase, iContrast} = buildColorGrating(pixelsPerDegree, [stimSize stimSize], ...
            spatialFreq, orientation, phase, contrast, ...
            1, 'bw', 1, 1);
        
        stim{iPhase, iContrast} = maskWithAnnulus(s{iPhase,iContrast}, ...
            length(s{iPhase,iContrast}), ...
            0, blurRadius, backgroundColor);
    end
end

% arrange images on the two sides of the screen in all possible
% combinations
% make spacer matrix
spacer = ones(size(stim{1,1},1), round(stimSpacerWidth*pixelsPerDegree)).*backgroundColor;

% generate all combinations of phase, contrast, and side of screen
for iP1 = 1:numel(phases)
    for iC1 = 1:numel(contrasts)
        for iP2 = 1:numel(phases)
            for iC2 = 1:numel(contrasts)
                
                stimMatrix{iP1,iC1,iP2,iC2} = ...
                    [stim{iP1,iC1} spacer stim{iP2,iC2}];
                
                stimIDs(iP1,iC1,iP2,iC2) = ...
                    iP1*1000 + iC1*100 + iP2*10 + iC2; % [phase left, contrast left, phase right, contrast right]
            end
        end
    end
end

% put into images, keeping track of the IDs
for iIm = 1:numel(stimMatrix)
    images(:,:,iIm) = stimMatrix{iIm};
    imageIDs(iIm,1) = stimIDs(iIm);
end
imageIDHeaders = {'left-phase', 'left-contrast', 'right-phase', 'right-contrast'};

% add blank stimulus
images(:,:,end+1) = ones(size(images(:,:,1)))*backgroundColor;
imageIDs(end+1) = 0;

% % show images (for debugging)
% for iIm = 1:size(images,3)
%     imshow(images(:,:,iIm));
%     pause(1);
% end

%% Generate the stimulus sequence
runDur = blockDur*nBlocks;
blockStartTimes = 0:blockDur:runDur-blockDur;
nFramesPerBlock = blockDur*refrate;

% randomly selected target times
targetStartTimes = [];
for iBlock = 1:nBlocks
    if ~strcmp(blockNames{blockOrder(iBlock)},'blank')
        nTargets = ceil(rand*maxTargetsPerBlock); % 1 to max targets
        targetTimes = choosePseudoRandomTargetTimes(blockDur-targetLeadTime-targetEndTime, nTargets, targetCushion);
        targetTimes = targetTimes + blockStartTimes(iBlock) + targetLeadTime;
        targetStartTimes = [targetStartTimes; targetTimes'];
    end
end
targetEndTimes = targetStartTimes + targetDur;
nTargets = numel(targetStartTimes);
% randomly assign targets to either the left or right stimulus
targetSides = mod(randperm(nTargets),2)+1;

% time points
seqtiming = (0:1/refrate:runDur-1/refrate)';

% phase sequences
fastPhaseSeq = repmat(fastUnit,1,ceil(nFramesPerBlock/numel(fastUnit)))';
slowPhaseSeq = repmat(slowUnit,1,ceil(nFramesPerBlock/numel(slowUnit)))';

% specify what's happening on every frame
blockIdx = 1;
blockName = blockNames{blockOrder(blockIdx)};
attBlockName = attBlockNames{attBlockOrder(blockIdx)};
phaseSeqIdx = 1;
targetIdx = 1;
targetOn = 0;
for iFrame = 1:numel(seqtiming)
    time = seqtiming(iFrame);
    % start a new block when it's time
    if blockIdx < nBlocks && ...
            time >= blockStartTimes(blockIdx+1) - 0.00001;
        blockIdx = blockIdx+1;
        blockName = blockNames{blockOrder(blockIdx)};
        attBlockName = attBlockNames{attBlockOrder(blockIdx)};
        phaseSeqIdx = 1;
        newBlock = 1; % used to decide whether block triggers will be on
    else
        newBlock = 0;
    end
    if iFrame==1 % make sure newBlock is 1 on the first frame
        newBlock = 1;
    end
    
    % determine if target is on
    % turn it on when it is time and leave on until time to turn off
    if targetIdx <= nTargets && targetStartTimes(targetIdx)-time < 1/refrate
        targetOn = 1;
    end
    if targetIdx <= nTargets && targetEndTimes(targetIdx)-time < 1/refrate
        targetOn = 0;
        targetIdx = targetIdx + 1;
    end
    if targetOn
        if targetSides(targetIdx)==1
            c1 = 2;
            c2 = 1;
            targetOnSeq(iFrame,1) = 1;
        elseif targetSides(targetIdx)==2
            c1 = 1;
            c2 = 2;
            targetOnSeq(iFrame,1) = 2;
        end
    else
        c1 = 1; c2 = 1;
        targetOnSeq(iFrame,1) = 0;
    end
    
    % determine image
    switch blockName
        case 'blank'
            p1 = 0; p2 = 0; c1 = 0; c2 = 0;
        case 'fast-left'
            p1 = fastPhaseSeq(phaseSeqIdx); % left
            p2 = slowPhaseSeq(phaseSeqIdx); % right
            %             c1 = 1; c2 = 1;
        case 'slow-left'
            p1 = slowPhaseSeq(phaseSeqIdx); % left
            p2 = fastPhaseSeq(phaseSeqIdx); % right
            %             c1 = 1; c2 = 1;
        otherwise
            error('blockName not recognized')
    end
    phaseSeqIdx = phaseSeqIdx + 1;
    imageID = 1000*p1 + 100*c1 + 10*p2 + c2;
    
    seq(iFrame,1) = find(imageIDs==imageID);
    
    % determine cue
    switch attBlockName
        case 'no-att'
            if time-blockStartTimes(blockIdx) < respDur
                % give a response window at the beginning of blank blocks
                fixSeq(iFrame,1) = 2;
            elseif blockIdx < nBlocks && (blockStartTimes(blockIdx+1)-time < attCueLeadTime)
                % cue the next attention block right before it starts
                switch attBlockNames{attBlockOrder(blockIdx+1)}
                    case 'att-left'
                        fixSeq(iFrame,1) = 4;
                    case 'att-right'
                        fixSeq(iFrame,1) = 5;
                end
            else
                fixSeq(iFrame,1) = 1;
            end
        case 'att-left'
            fixSeq(iFrame,1) = 4;
        case 'att-right'
            fixSeq(iFrame,1) = 5;
        otherwise
            error('attBlockName not recognized')
    end
    
    switch triggerOption
        case 'conditionID'
            % determine trigger - condition ID
            if newBlock % only give condition trigger at the first frame of the block
                switch blockName
                    case 'blank'
                        trig = 7;
                    case 'fast-left'
                        if strcmp(attBlockName,'att-left')
                            trig = 1;
                        elseif strcmp(attBlockName,'att-right')
                            trig = 2;
                        end
                    case 'slow-left'
                        if strcmp(attBlockName,'att-left')
                            trig = 3;
                        elseif strcmp(attBlockName,'att-right')
                            trig = 4;
                        end
                    otherwise
                        error('blockName not recognized')
                end
%                 blockTrig = trig;
            elseif targetOnSeq(iFrame)~=0
                % triger for target side, only on first target frame
                if targetOnSeq(iFrame)==1 && targetOnSeq(iFrame-1)==0
                    trig = 5; % target on left
                elseif targetOnSeq(iFrame)==2 && targetOnSeq(iFrame-1)==0
                    trig = 6; % target on right
                else
                    trig = NaN;
                end
%                 targetTrig = trig;
            else
                trig = NaN;
%                 blockTrig = NaN;
%                 targetTrig = NaN;
            end
            
            trigSeq(iFrame,1) = computeTrigger(trig);
%             trigSeq(iFrame,1) = computeTrigger(blockTrig, targetTrig);
%             trigSeq(iFrame,1) = 0;
            
        case 'combinatorial'
            % determine trigger - combinatorial code
            if newBlock % only give condition trigger at the first frame of the block
                if strcmp(blockName,'blank')
                    blankTrig = 7;
                else
                    blankTrig = NaN;
                end
                if strcmp(blockName,'fast-left')
                    fastSideTrig = 1;
                elseif strcmp(blockName,'slow-left')
                    fastSideTrig = 2;
                else
                    fastSideTrig = NaN;
                end
                if strcmp(attBlockName,'att-left')
                    attSideTrig = 3;
                elseif strcmp(attBlockName,'att-right')
                    attSideTrig = 4;
                else
                    attSideTrig = NaN;
                end
            else
                blankTrig = NaN;
                fastSideTrig = NaN;
                attSideTrig = NaN;
            end
            
            % triger for target side, only on first target frame
            if targetOnSeq(iFrame)==1 && targetOnSeq(iFrame-1)==0
                targetTrig = 5; % target on left
            elseif targetOnSeq(iFrame)==2 && targetOnSeq(iFrame-1)==0
                targetTrig = 6; % target on right
            else
                targetTrig = NaN;
            end
            
            trigSeq(iFrame,1) = computeTrigger(blankTrig, fastSideTrig, attSideTrig, targetTrig);
%             trigSeq(iFrame,1) = 0;
 
        otherwise
            error('triggerOption not recognized')
    end
end

% show targetOnSeq and trigSeq
f(1) = figure;
subplot(2,1,1)
plot(seqtiming,targetOnSeq)
subplot(2,1,2)
plot(seqtiming,trigSeq)

% display triggers by channel
f(2) = displayTrigger(trigSeq, nBlocks);
set(f(2),'Position',[0 0 1200 900]);

%% Create stimulus strucutre
% set remaining stimulus variables
cmap = repmat((0:255)',1,3);
srcRect = [0 0 size(images,2) size(images,1)];
destRect = CenterRectOnPoint(srcRect, cx, cy+stimPos(2)*pixelsPerDegree);
diodeSeq = repmat([0 0 1 1], 1, ceil(length(seq)/4))';

% store in stimulus structure
stimulus.images = images*255;
stimulus.cmap = cmap;
stimulus.seq = seq;
stimulus.seqtiming = seqtiming;
stimulus.fixSeq = fixSeq;
stimulus.srcRect = srcRect;
stimulus.destRect = destRect;
stimulus.trigSeq = trigSeq;
stimulus.diodeSeq = diodeSeq;

% save stimulus
if saveStim
    save(sprintf('%s/%s.mat', stimDir, stimFile), 'stimulus', 'p')
end

% save figs
if saveFigs
    rd_saveAllFigs(f, {'trigplot','trigchan'}, stimFile);
end

