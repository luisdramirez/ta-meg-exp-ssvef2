function makeTADetectDiscrimStim(run)

%% run setup
% run = 3;
saveStim = 0;
saveFigs = 0;

%% add paths
addpath(genpath('/../vistadisp'))
addpath('../TAPilot')


%% file i/o
stimDir = 'stimuli';
stimFile = sprintf('taDetect%d', run);

%% screen setup
displayName = 'meg_lcd';
d = loadDisplayParams(displayName);
pixelsPerDegree = 1/d.pixelSize;
screenWidth = d.numPixels(1); % (px)
screenHeight = d.numPixels(2); % (px)
cx = round(screenWidth/2);
cy = round(screenHeight/2);

%% keys setup
keyNames = {'1!','2@'}; % [present absent]
keyCodes = KbName(keyNames);

%% timing setup
refrate = 60; % (Hz)
blockDur = 5; % (s)
targetDur = 3/refrate; % (s)
targetLeadTime = 1.5; % (s) % no targets in first part of block
targetSOA = 0.8; % (s) % SOA between targets
cueTargetSOA = 1; % (s) % SOA between cues and targets, same for pre- and post-cues
attCueLeadTime = 0.5; % (s)
respDur = 1.4; % (s)
feedbackDur = 0.3;
cueDur = 0.1;
if refrate==75
    % 75 Hz SSVEP unit sequences: 4 frames (75/4=18.75 Hz) and 5 frames (75/5=15 Hz)
    fastUnit = [1 1 2 2]; % gives the phase (1 or 2) of each frame
    slowUnit = [1 1 1 2 2];
else
    % 60 Hz SSVEP unit sequences: 3 frames (60/3=20 Hz) and 4 frames (60/4=15 Hz)
    fastUnit = [1 2 2]; % gives the phase (1 or 2) of each frame
    slowUnit = [1 1 2 2];
end

%% blocks setup (one run)
blockNames = {'blank','fast-left'}; % fast-left
attBlockNames = {'no-att','att-right'}; % att-right
targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres','abs-abs'};
cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'}; % 2-1 = cueT2,postcueT1
[blockOrder,attBlockOrder, targetBlockOrder,cueBlockOrder] = block_gen(blockNames,attBlockNames, targetBlockNames ,cueBlockNames );
nBlocks = numel(blockOrder);

%% stim setup   
stimSize = 8;
spatialFreq = 1;
orientation = 0;
stimContrast = 0.64;
targetContrast = 0.84;
contrasts = [stimContrast targetContrast];
blurRadius = 0.2;
backgroundColor = 128/255;
phases = [0 pi];

stimPos = [6 4]; % [x y]
stimSpacerWidth = (stimPos(1)-stimSize/2)*2;

%% sound setup
Fs = 44100;
cueFreqs = [784 523];
% cueFreqs = [1046.5 440]; % [higher high C = target 1, lower A = target 2]
for iTone = 1:numel(cueFreqs)
    tone = MakeBeep(cueFreqs(iTone), cueDur, Fs);
    cueTones(iTone,:) = applyEnvelope(tone, Fs);
end

%% trigger setup
triggerOption = 'conditionID'; % 'conditionID','combinatorial'

%% Store all stim params
% sorry this is kind of a horrible way to do this
p = v2struct(...
    displayName, pixelsPerDegree, screenWidth, screenHeight, cx, cy, ...
    keyNames, keyCodes, ...
    refrate, blockDur, targetDur, targetLeadTime, targetSOA, ...
    attCueLeadTime, respDur, feedbackDur, fastUnit, slowUnit, ...
    blockNames, blockOrder, attBlockNames, attBlockOrder, targetBlockNames, targetBlockOrder, ...
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

%% Determine the stimulus times
runDur = blockDur*nBlocks;
blockStartTimes = 0:blockDur:runDur-blockDur;
nFramesPerBlock = blockDur*refrate;

% fixed target times (T1 and T2) on the attended (right-side) stimulus
targetStartTimes = [];
targetAbsStartTimes = [];
cueStartTimes = [];
for iBlock = 1:nBlocks
    if ~strcmp(blockNames{blockOrder(iBlock)},'blank')
        targetTimes = [0 targetSOA];
        cueTimes = [targetTimes(1)-cueTargetSOA targetTimes(2)+cueTargetSOA];
        
        % eliminate "absent" targets, depending on condition
        targetBlockName = targetBlockNames{targetBlockOrder(iBlock)};
        switch targetBlockName
            case 'pres-pres'
                targetAbsTimes = [];
            case 'pres-abs'
                targetAbsTimes = targetTimes(2);
                targetTimes(2) = [];
            case 'abs-pres'
                targetAbsTimes = targetTimes(1);
                targetTimes(1) = [];
            case 'abs-abs'
                targetAbsTimes = targetTimes;
                targetTimes = [];
            otherwise
                error('targetBlockName not recognized')
        end
        
        % target present
        targetTimes = targetTimes + blockStartTimes(iBlock) + targetLeadTime;
        targetStartTimes = [targetStartTimes; targetTimes'];
        
        % target absent
        targetAbsTimes = targetAbsTimes + blockStartTimes(iBlock) + targetLeadTime;
        targetAbsStartTimes = [targetAbsStartTimes; targetAbsTimes'];
        
        % pre- and post-cues
        cueTimes = cueTimes + blockStartTimes(iBlock) + targetLeadTime;
        cueStartTimes = [cueStartTimes; cueTimes'];
    end
end
targetEndTimes = targetStartTimes + targetDur;
nTargets = numel(targetStartTimes);
% all these targets are on the right
targetSides = ones(1,nTargets)*2;

% target absent is treated in the same way
targetAbsEndTimes = targetAbsStartTimes + targetDur;
nAbsTargets = numel(targetAbsStartTimes);
targetAbsSides = ones(1,nAbsTargets)*2;

% tone cues
nCues = numel(cueStartTimes);

%% Generate the stimulus sequence
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
targetIdx = 1; % target present (on either side)
targetOn = 0;
targetAbsIdx = 1; % target absent
targetAbsOn = 0;
cueIdx = 1;
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
    
    % determine if is time to turn on the tone cue
    if cueIdx <= nCues && cueStartTimes(cueIdx)-time < 1/refrate
        cueBlockName = cueBlockNames{cueBlockOrder(blockIdx)};
        if mod(cueIdx,2) % odd cues are pre-cues
            cueType = str2double(cueBlockName(1));
        else % even cues are response cues
            cueType = str2double(cueBlockName(end));
        end
        if isnan(cueType)
            error('cues should not be presented during no-cue blocks. check code.')
        end
        cueSeq(iFrame,1) = cueType;
        cueIdx = cueIdx + 1;
    else
        cueSeq(iFrame,1) = 0;
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
    
    % determine if an absent target is "on" (if it is time for that target)
    % treat it the same way as real targets, in order to determine
    % triggers, but do not change the stimulus contrast
    if targetAbsIdx <= nAbsTargets && targetAbsStartTimes(targetAbsIdx)-time < 1/refrate
        targetAbsOn = 1;
    end
    if targetAbsIdx <= nAbsTargets && targetAbsEndTimes(targetAbsIdx)-time < 1/refrate
        targetAbsOn = 0;
        targetAbsIdx = targetAbsIdx + 1;
    end
    if targetAbsOn
        if targetAbsSides(targetAbsIdx)==1
            targetAbsOnSeq(iFrame,1) = 1;
        elseif targetSides(targetAbsIdx)==2
            targetAbsOnSeq(iFrame,1) = 2;
        end
    else
        targetAbsOnSeq(iFrame,1) = 0;
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
    
    % determine spatial attention cue
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
            if blockStartTimes(blockIdx+1)-time < feedbackDur
                % display feedback at the end of the block
                fixSeq(iFrame,1) = 8; % blue
             else
                fixSeq(iFrame,1) = 4;
            end
        case 'att-right'
            if blockStartTimes(blockIdx+1)-time < feedbackDur
                % display feedback at the end of the block
                fixSeq(iFrame,1) = 8; % blue
            else
                fixSeq(iFrame,1) = 5;
            end
        otherwise
            error('attBlockName not recognized')
    end
    
    % determine correct response keycode
    cueBlock = cueBlockNames{cueBlockOrder(blockIdx)};
    if ~strcmp(cueBlock, 'no-cue') && ...
            time-blockStartTimes(blockIdx) > targetLeadTime + targetSOA + cueTargetSOA && ...
            blockStartTimes(blockIdx+1)-time > feedbackDur
        % which target is post-cued?
        responseCue = str2double(cueBlock(end));
        switch targetBlockNames{targetBlockOrder(blockIdx)}
            case 'pres-pres'
                correctResponse = 1; % present
            case 'pres-abs'
                if responseCue==1
                    correctResponse = 1; % present
                else
                    correctResponse = 2; % absent
                end
            case 'abs-pres'
                if responseCue==1
                    correctResponse = 2; % absent
                else
                    correctResponse = 1; % present
                end
            case 'abs-abs'
                correctResponse = 2; % absent
        end
        keyCodeSeq(iFrame,1) = keyCodes(correctResponse);
    else
        keyCodeSeq(iFrame,1) = 0;
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
                            trig = NaN; % never happens, so don't use up the trigger
                        elseif strcmp(attBlockName,'att-right')
                            trig = 1;
                        end
                    case 'slow-left'
                        if strcmp(attBlockName,'att-left')
                            trig = NaN; % never happens, so don't use up the trigger
                        elseif strcmp(attBlockName,'att-right')
                            trig = 2;
                        end
                    otherwise
                        error('blockName not recognized')
                end
            elseif targetOnSeq(iFrame)~=0
                % triger for target side, only on first target frame
                if targetOnSeq(iFrame)==1 && targetOnSeq(iFrame-1)==0
                    trig = 5; % target on left
                elseif targetOnSeq(iFrame)==2 && targetOnSeq(iFrame-1)==0
                    trig = 6; % target on right
                else
                    trig = NaN;
                end
            elseif targetAbsOnSeq(iFrame)~=0
                % triger for target side, only on first target frame
                if targetAbsOnSeq(iFrame)==1 && targetAbsOnSeq(iFrame-1)==0
                    trig = 3; % target on left
                elseif targetAbsOnSeq(iFrame)==2 && targetAbsOnSeq(iFrame-1)==0
                    trig = 4; % target on right
                else
                    trig = NaN;
                end
            elseif cueSeq(iFrame)~=0
                trig = 8; % pre- or post-cue
            else
                trig = NaN;
            end
            
            trigSeq(iFrame,1) = computeTrigger(trig);
            
        case 'combinatorial'
            error('sorry, combinatorial trigger code not yet implemented')
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
subplot(3,1,1)
plot(seqtiming,targetOnSeq)
subplot(3,1,2)
plot(seqtiming,trigSeq)
subplot(3,1,3)
plot(seqtiming,keyCodeSeq)

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
stimulus.sounds = cueTones;
stimulus.cmap = cmap;
stimulus.seq = seq;
stimulus.seqtiming = seqtiming;
stimulus.fixSeq = fixSeq;
stimulus.srcRect = srcRect;
stimulus.destRect = destRect;
stimulus.trigSeq = trigSeq;
stimulus.diodeSeq = diodeSeq;
stimulus.keyCodeSeq = keyCodeSeq;
stimulus.soundSeq = cueSeq;

% store in order structure
order.blockorder = blockOrder;
order.attBlockOrder = attBlockOrder;
order.targetBlockOrder = targetBlockOrder;
order.cueBlockOrder = cueBlockOrder;

% save stimulus
if saveStim
    save(sprintf('%s/%s.mat', stimDir, stimFile), 'stimulus', 'p','order')
end

% save figs
if saveFigs
    rd_saveAllFigs(f, {'trigplot','trigchan'}, stimFile);
end