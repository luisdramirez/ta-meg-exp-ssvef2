function makeTADetectDiscrimStim(subjectID, run, displayName)

%% run setup
%run = 7;
saveStim = 1;
saveFigs = 0;   

%% add paths
addpath(genpath('../../vistadisp'))
addpath('../TAPilot')

%% file i/o
stimDir = 'stimuli';
vistaStimPath = 'vistadisp/Applications2/Retinotopy/standard/storedImagesMatrices';
vistaStimDir = sprintf('../../%s', vistaStimPath);
stimFile = sprintf('%s_taDetectDiscrim%d', subjectID, run);

%% screen setup
if nargin < 3
    displayName = 'meg_lcd'; % 'meg_lcd','Carrasco_L2','Carrasco_L1','Carrasco_R1'
end
d = loadDisplayParams(displayName);

% pixelsPerDegree = 1/d.pixelSize;
pixelsPerCm = 1/d.pixelSize;
degPerCm = atan(1/d.distance)*180/pi;
pixelsPerDegree = pixelsPerCm/degPerCm;

screenWidth = d.numPixels(1); % (px)
screenHeight = d.numPixels(2); % (px)
cx = round(screenWidth/2);
cy = round(screenHeight/2);

%% keys setup
responseOption = 'targetContrast4Levels'; % 'targetType','targetPos','targetContrast4Levels'
% keyNames = {'2@','3#'}; % [contrast2 contrast3]
keyNames = {'1!','2@','3#','4$'}; % [target1 target2 absent] or [contrast1 contrast2 contrast3 contrast4]
keyCodes = KbName(keyNames);

%% timing setup
refrate = 60; % (Hz)
nFramesPerTarget = 6; % 8
targetDur = nFramesPerTarget/refrate; % (s)
targetLeadTime = 1.5; % (s) % no targets in first part of block
targetSOA = 18/60; %18/60 = .300; 15/60 = .250; 16/60 = .267; %0.6; % (s) % SOA between targets (- difference from .8)
cueTargetSOA = 1; % (s) % SOA between cues and targets, same for pre- and post-cues
attCueLeadTime = 0.5; % (s)
respDur = 1.6; %1.2; % (s) % if unlimited response window, then 1 frame (see showScanStimulus)
feedbackDur = 0.3; % (s)
cueDur = 0.1; % (s)
blockDur = targetLeadTime + targetSOA + cueTargetSOA + respDur + feedbackDur; % (s)
iti = 1.5; %2; % (s)
jitter = 'blockPrecueInterval'; % 'blockPrecueInterval', 'ITI', 'none' % add jittered interval between trials
flickerType = 'counterphase'; % 'counterphase','onoff'
if refrate==60
    % 60 Hz SSVEP unit sequences: 3 frames (60/3=20 Hz) and 4 frames (60/4=15 Hz)
    switch flickerType
        case 'counterphase'
            fastUnit = [1 2 2]; % gives the phase (1 or 2) of each frame
            slowUnit = [1 1 1 2 2 2]; % [1 1 2 2] = 30 Hz;
        case 'onoff'
            fastUnit = [1 1 0];
            slowUnit = [1 0 1 0];
        otherwise
            error('flickerType not recognized')
    end
else
    error('refrate must be 60 Hz')
end

%% target setup
target.type = 'contrast'; % 'dot','lines','grating','cb','contrast'
target.catchTrials = false;

%% blocks setup (one run)
blockNames = {'blank','fast-left'}; % fast-left, slow-left
attBlockNames = {'no-att','att-right'}; % att-right
targetBlockNames = {'no-targ','pres-pres'};
% targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres','abs-abs'};
cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'}; % 2-1 = cueT2,postcueT1
% cueBlockNames = {'no-cue','1-1','2-2'};
% [blockOrder, attBlockOrder, targetBlockOrder, cueBlockOrder, targetTypeBlockOrder] ...
%     = block_gen(blockNames,attBlockNames, targetBlockNames, cueBlockNames, run, target.catchTrials);
[blockOrder, attBlockOrder, targetBlockOrder, cueBlockOrder, targetTypeBlockOrder, targetPedestalBlockOrder] = ...
    block_gen3(subjectID, run, stimDir);

%%% debugging %%%
% blockOrder = blockOrder(1:6);
% attBlockOrder = attBlockOrder(1:6);
% targetBlockOrder = targetBlockOrder(1:6);
% cueBlockOrder = cueBlockOrder(1:6);
% targetTypeBlockOrder = targetTypeBlockOrder(:,1:6);
%%%%%

nBlocks = numel(blockOrder);

if target.catchTrials
    targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres'};
end

%% stim setup  
stimType = 'radialcb'; %'grating' 'checkerboard' 'bullseye' 'radialcb' 'spiralcb' 'radialcbgrad'
stimSize = 4;
spatialFreq = 3;
orientation = 0;
possibleContrasts = [
    0.2812
    0.2922
    0.3036
    0.3155
    0.3278
    0.3407
    0.3540
    0.3678
    0.3822
    0.3972
    0.4127
    0.4289
    0.4457
    0.4631
    0.4812
    0.5000
    0.5196
    0.5399
    0.5610
    0.5830
    0.6058
    0.6295
    0.6541
    0.6797
    0.7063
    0.7340
    0.7627
    0.7925
    0.8235
    0.8557
    0.8892];
stimContrast = 0.4; % 0.64
targetContrast = 0.4; % 0.64
contrasts = [stimContrast targetContrast];
blurRadius = 0.2;
backgroundColor = 128/255;
phases = [0 pi];

radialCB.thetaCycles = 8;
radialCB.A = 1;
switch stimType
    case {'radialcb','radialcbgrad'}
        radialCB.b = 0.2;
        radialCB.E = 0.05;
    case 'spiralcb'
        radialCB.b = 0.4;
        radialCB.E = 0.1;
end
radialCB.gradientAngles = [-135 -45 135 45]; % top left higher contrast, top right, bottom left, bottom right

% fixation
fixDiam = 0.15;

stimPos = [0 0]; % [x y]
stimSpacerWidth = (stimPos(1)-stimSize/2)*2;

%% sound setup
Fs = 44100;
% cueFreqs = [784 523];
% cueFreqs = [1300 250] % antonio
cueFreqs = [1046.5 440]; % [higher high C = target 1, lower A = target 2]
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
    responseOption, keyNames, keyCodes, ...
    refrate, blockDur, targetDur, targetLeadTime, targetSOA, cueTargetSOA, ...
    attCueLeadTime, respDur, feedbackDur, fastUnit, slowUnit, ...
    blockNames, blockOrder, attBlockNames, attBlockOrder, targetBlockNames, targetBlockOrder, ...
    cueBlockNames, cueBlockOrder, targetTypeBlockOrder, ...
    stimSize, stimPos, spatialFreq, orientation, stimContrast, targetContrast, ...
    contrasts, blurRadius, backgroundColor, phases, radialCB, triggerOption, jitter, flickerType);

%% Make the stimuli
for iPhase = 1:numel(phases)
    for iContrast = 1:numel(contrasts)
        phase = phases(iPhase);
        contrast = contrasts(iContrast);
        
        switch stimType
            case 'grating'
                s{iPhase, iContrast} = buildColorGrating(pixelsPerDegree, [stimSize stimSize], ...
                    spatialFreq, orientation, phase, contrast, ...
                    1, 'bw', 1, 1);
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
                s{iPhase, iContrast} = (c-0.5)*contrast+0.5;
            case 'bullseye'
                s{iPhase, iContrast} = CreateSpiral(d, stimSize, spatialFreq, phase, contrast)./2 + .5;
                %s{iPhase, iContrast} = (c-0.5)*contrast+0.5;
            case {'radialcb','radialcbgrad'}
                s{iPhase, iContrast} = makeRadialCheckerboard(pixelsPerDegree, stimSize, phase, contrast, ...
                    radialCB.thetaCycles, radialCB.E, radialCB.A, radialCB.b);
            case 'spiralcb'
                s{iPhase, iContrast} = makeSpiralCheckerboard(pixelsPerDegree, stimSize, phase, contrast, ...
                    radialCB.thetaCycles, radialCB.E, radialCB.A, radialCB.b);
            otherwise
                error('stimType not recognized')
        end
        
        stim{iPhase, iContrast} = maskWithAnnulus(s{iPhase,iContrast}, ...
            length(s{iPhase,iContrast}), ...
            0, blurRadius, backgroundColor); %stimulus image generated here
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
                
%                 stimMatrix{iP1,iC1,iP2,iC2} = ...
%                     [stim{iP1,iC1} spacer stim{iP2,iC2}];
                stimMatrix{iP2,iC2} = ... %iP1,iC1,
                    [stim{iP2,iC2}]; %stim{iP1,iC1} 
                
%                 stimIDs(iP1,iC1,iP2,iC2) = ...
%                     iP1*1000 + iC1*100 + iP2*10 + iC2; % [phase left, contrast left, phase right, contrast right]
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

% % show images (for debugging)
% for iIm = 1:size(images,3)
%     imshow(images(:,:,iIm));
%     pause(1);
% end

%% Set up targets
% set up lines, just vertical and horizontal of the right size, centered on
% (0,0)
xy0 = round([0 0 stimSize/2 -stimSize/2; ...
    stimSize/2 -stimSize/2 0 0].*pixelsPerDegree*0.85); % so it doesn't go all the way to the edge of the patch
targetCenter = [cx cy] + stimPos.*pixelsPerDegree;

% store in target structure
switch target.type
    case 'lines'
        target.center = targetCenter;
        target.colors = [1 0 0]*255; % red
        target.xy0 = xy0;
        target.width = 2;
        target.baseOrient = 45;
    case 'dot'
        target.center = targetCenter;
        target.colors = [1 1 1]*255/2; % gray
        target.maxRadiusPx = stimSize/2*pixelsPerDegree*0.85;
        target.pixelsPerDegree = pixelsPerDegree;
        target.dotLocs = [1 -1]; % lower or upper part of target
    case 'grating'
        target.contrast = contrasts(1);
        target.phases = phases;
        target.pixelsPerDegree = pixelsPerDegree;
        target.stimSize = stimSize;
        target.spatialFreq = spatialFreq;
        target.orientation = orientation;
        target.blurRadius = blurRadius;
        target.backgroundColor = backgroundColor;
        target.spacer = spacer;
        target.stim = stim;
    case 'cb'
        target.pixelsPerDegree = pixelsPerDegree;
        target.imSize = stimSize; % whole grating square
        target.stimSize = stimSize;
        target.size = 1.5; % 0.5 % sigma of gaussian aperture
        target.spatialFreq = 4;
        target.center = targetCenter;
    case 'contrast'
        % SPECIFY TARGET PARAMETERS HERE
        target.phases = phases;
        target.stimType = stimType;
        target.cx = cx;
        target.cy = cy;
        target.pixelsPerDegree = pixelsPerDegree;
        target.stimSize = stimSize;
        target.spatialFreq = spatialFreq;
        target.orientation = orientation;
        target.blurRadius = blurRadius;
        target.backgroundColor = backgroundColor;
        target.backgroundContrast = stimContrast;
        target.stim = stim;
        target.radialCB = radialCB;
        target.nFramesPerTarget = nFramesPerTarget;
        target.positions = 1; % (1:8)';
        target.nPedestals = 2; % 1 for normal discrimination from baseline
        if numel(target.positions)==1
            target.sigma = target.stimSize*2;
        else
            target.sigma = target.stimSize/8;
        end
        % CREATE TARGET POSITIONS 
        nTargetAppears = length(targetBlockOrder(targetBlockOrder == 2))*2 + ...
            length(targetBlockOrder(targetBlockOrder == 3)) +...
            length(targetBlockOrder(targetBlockOrder == 4));
        if any(strcmp(targetBlockNames,'pres-abs')) && ~target.catchTrials
            positions_mat = repmat(target.positions, 1, nTargetAppears/length(target.positions)); 
            
            % Pedestals (for contrast discrimination)
            pedestals_mat = repmat((1:target.nPedestals)', 1, nTargetAppears/target.nPedestals);
        else
            if target.catchTrials
                positions_mat = repmat(target.positions, nTargetAppears/length(target.positions)/2, 4);
                pedestals_mat = repmat((1:target.nPedestals)', nTargetAppears/target.nPedestals/2, 4);
            else
                positions_mat = repmat(target.positions, nTargetAppears/length(target.positions)/2, 2);
                pedestals_mat = repmat((1:target.nPedestals)', nTargetAppears/target.nPedestals/2, 2);
            end
        end
        posShuffled = Shuffle(positions_mat);
        posShuffledHeaders = {'pres-presT1', 'pres-presT2','pres-abs','abs-pres'};
        pedestalShuffled = Shuffle(pedestals_mat);
        
        % Generate guassian center coordinates 
        if numel(target.positions)==1
            target.coords = [0 0];
        else
            xmax = size(target.stim{1},1); ymax = size(target.stim{1},2);
            cx2 = 0; cy2 = 0;  %origin of coordiantes
            r = xmax/4; %radius of center coordinates
            % start from 6 o'clock and go clockwise (left then right)
            theta = 22.5+90:45:360+90; theta = deg2rad(mod(theta,360));
            x0 = cx2 + r * cos(theta); %generate x coordinates
            y0 = cx2 + r * sin(theta); %generate y coordinates
            gaussCoords = [x0' y0']; %store coordinates
            target.coords = gaussCoords; %store coordinates
            target.responsePosSets = ([1:numel(target.positions)/2; numel(target.positions)/2+1:numel(target.positions)])';
        end
        
    otherwise
        error('target.type not recognized')
end

%% Determine the stimulus times
switch jitter
    case 'blockPrecueInterval'
        itiSeq = ones(1,nBlocks)*iti;
        blockDur = blockDur + iti;
        jit = (0:0.2:1); % recall there is always 0.5 s before cue
        jitSeq = Shuffle(repmat(jit,1,ceil(nBlocks/numel(jit))));
        jitSeq = jitSeq(1:nBlocks);
        runDur = blockDur*nBlocks + sum(jitSeq);
        blockStartTimes = (0:blockDur:blockDur*nBlocks-blockDur) + cumsum([0 jitSeq(1:end-1)]);
        nFramesPerBlock = (blockDur + max(jitSeq))*refrate; % number of frames in the longest block
    case 'ITI'
        iti = 0:0.2:1; % recall there is always 0.5 s before cue
        itiSeq = Shuffle(repmat(iti,1,ceil(nBlocks/numel(iti))));
        itiSeq = itiSeq(1:nBlocks);
        runDur = blockDur*nBlocks + sum(itiSeq);
        blockStartTimes = (0:blockDur:blockDur*nBlocks-blockDur) + cumsum([0 itiSeq(1:end-1)]);
        nFramesPerBlock = (blockDur + max(itiSeq))*refrate; % number of frames in the longest block
    case 'none'
        itiSeq = ones(1,nBlocks)*iti;
        blockDur = blockDur + iti;
        runDur = blockDur*nBlocks;
        blockStartTimes = 0:blockDur:runDur-blockDur;
        nFramesPerBlock = blockDur*refrate;
    otherwise
        error('jitter not recognized')
end
p.blockDur = blockDur; % as this may have changed

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
        targetTimes = targetTimes + blockStartTimes(iBlock) + targetLeadTime + jitSeq(iBlock);
        targetStartTimes = [targetStartTimes; targetTimes'];
        
        % target absent
        targetAbsTimes = targetAbsTimes + blockStartTimes(iBlock) + targetLeadTime + jitSeq(iBlock);
        targetAbsStartTimes = [targetAbsStartTimes; targetAbsTimes'];
        
        % pre- and post-cues
        cueTimes = cueTimes + blockStartTimes(iBlock) + targetLeadTime + jitSeq(iBlock);
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
posRowCount = ones(1,size(posShuffled,2));

for iFrame = 1:numel(seqtiming)
    time = seqtiming(iFrame);
    % start a new block when it's time
    if blockIdx < nBlocks && ...
            time >= blockStartTimes(blockIdx+1) - 0.00001;
        blockIdx = blockIdx+1;
        blockName = blockNames{blockOrder(blockIdx)};
        attBlockName = attBlockNames{attBlockOrder(blockIdx)};
        cueBlockName = cueBlockNames{cueBlockOrder(blockIdx)};
        targetBlockName = targetBlockNames{targetBlockOrder(blockIdx)};
        phaseSeqIdx = 1;
        newBlock = 1; % used to decide whether block triggers will be on
    else
        newBlock = 0;
    end
    if iFrame==1 % make sure newBlock is 1 on the first frame
        newBlock = 1;
    end
    
    % determine if is time to turn on the tone cue
    if cueIdx <= nCues && cueStartTimes(cueIdx)-time < 1/refrate - 0.00001
        cueBlockName = cueBlockNames{cueBlockOrder(blockIdx)};
        if mod(cueIdx,2) % odd cues are pre-cues
            cueType = str2double(cueBlockName(1));
            prePostCueSeq(iFrame,1) = 1;
        else % even cues are response cues
            cueType = str2double(cueBlockName(end));
            prePostCueSeq(iFrame,1) = 2;
        end
        if isnan(cueType)
            error('cues should not be presented during no-cue blocks. check code.')
        end
        cueSeq(iFrame,1) = cueType;
        cueIdx = cueIdx + 1;
    else
        cueSeq(iFrame,1) = 0;
        prePostCueSeq(iFrame,1) = 0;
    end
    
    % determine if target is on
    % turn it on when it is time and leave on until time to turn off
    if targetIdx <= nTargets && targetStartTimes(targetIdx)-time < 1/refrate - 0.00001
        targetOn = 1;
    end
    if targetIdx <= nTargets && targetEndTimes(targetIdx)-time < 1/refrate - 0.00001
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
    % specify the target type (for discrimination)
    if targetOn
        if targetOnSeq(iFrame-1)==0 % first frame of target
%             targetType = randi(2); % 1 or 2
            timeSinceCue = time - cueStartTimes(cueIdx-1);
            if timeSinceCue - (targetLeadTime-attCueLeadTime) < 1e-4
                whichTarget = 1;
            elseif timeSinceCue - (targetLeadTime-attCueLeadTime+targetSOA) < 1e-4
                whichTarget = 2;
            else
                error('target onset time does not match expected values')
            end
            targetType = targetTypeBlockOrder(whichTarget, blockIdx);
            targetTypeSeq(iFrame,1) = targetType;
            targetTypes(targetIdx) = targetType;
            
            targetPedestal = targetPedestalBlockOrder(whichTarget, blockIdx);
            pedestalSeq(iFrame,1) = targetPedestal;
            targetPedestals(targetIdx) = targetPedestal;
        else
            targetTypeSeq(iFrame,1) = targetTypeSeq(iFrame-1);
            pedestalSeq(iFrame,1) = pedestalSeq(iFrame-1);
        end
        switch targetBlockName
            case 'pres-pres'
                posCol = whichTarget;
            case 'pres-abs'
                posCol = 3;
            case 'abs-pres'
                posCol = 4;
            otherwise
                error('targetBlockName not viable')
        end
        
        posRow = posRowCount(posCol);
        posSeq(iFrame,1) = posShuffled(posRow,posCol);
%         pedestalSeq(iFrame,1) = pedestalShuffled(posRow,posCol);

%         if posSeq(iFrame-nFramesPerTarget+1)>0 %%% rd check w/ 9 frames
        if posSeq(iFrame-nFramesPerTarget+1)>0 && posRowCount(posCol) <= size(posShuffled,1)
            posRowCount(posCol) = posRowCount(posCol)+1;
        end
    else
        targetTypeSeq(iFrame,1) = 0;
        posSeq(iFrame,1) = 0;
        pedestalSeq(iFrame,1) = 0;
    end
    
    % determine if an absent target is "on" (if it is time for that target)
    % treat it the same way as real targets, in order to determine
    % triggers, but do not change the stimulus contrast
    if targetAbsIdx <= nAbsTargets && targetAbsStartTimes(targetAbsIdx)-time < 1/refrate - 0.00001
        targetAbsOn = 1;
    end
    if targetAbsIdx <= nAbsTargets && targetAbsEndTimes(targetAbsIdx)-time < 1/refrate - 0.00001
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
            % if in iti
            if blockIdx < nBlocks && (blockStartTimes(blockIdx+1)-time < itiSeq(blockIdx) - 0.00001)
                p1 = 0; p2 = 0;  c1 = 0; c2 = 0;
            else
                p1 = fastPhaseSeq(phaseSeqIdx); % left
                p2 = slowPhaseSeq(phaseSeqIdx); % right
                %             c1 = 1; c2 = 1;
            end
        case 'slow-left'
            if blockIdx < nBlocks && (blockStartTimes(blockIdx+1)-time < itiSeq(blockIdx) - 0.00001)
                p1 = 0; p2 = 0;  c1 = 0; c2 = 0;
            else
                p1 = slowPhaseSeq(phaseSeqIdx); % left
                p2 = fastPhaseSeq(phaseSeqIdx); % right
                %             c1 = 1; c2 = 1;
            end
        otherwise
            error('blockName not recognized')
    end
    if p1==0, c1 = 0; end
    if p2==0, c2 = 0; end
    phaseSeqIdx = phaseSeqIdx + 1;
    imageID = 10*p2 + c2; %(1000*p1 + 100*c1 + ) removed from front
    
    seq(iFrame,1) = find(imageIDs==imageID);
    
    if strcmp(target.type,'grating')
        target.phSeq(iFrame,:) = [p1 p2]; % keep track of phase of right grating
    end
    
    % determine spatial attention cue
    switch attBlockName
        case 'no-att'
%             if time-blockStartTimes(blockIdx) < respDur - 0.00001
%                 % give a response window at the beginning of blank blocks
%                 fixSeq(iFrame,1) = 2;
            if blockIdx < nBlocks && (blockStartTimes(blockIdx+1)-time < attCueLeadTime - 0.00001)
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
            if blockStartTimes(blockIdx+1)-time < feedbackDur + itiSeq(blockIdx) - 0.00001 && ...
                blockStartTimes(blockIdx+1)-time > itiSeq(blockIdx) - 0.00001
                % display feedback at the end of the block
                if target.catchTrials && ~isempty(strfind(targetBlockName,'abs'))
                    fixSeq(iFrame,1) = 9; % gray
                else
                    fixSeq(iFrame,1) = 8; % blue
                end
             else
                fixSeq(iFrame,1) = 4;
            end
        case 'att-right'
            if blockStartTimes(blockIdx+1)-time < feedbackDur + itiSeq(blockIdx) - 0.00001 && ...
                blockStartTimes(blockIdx+1)-time > itiSeq(blockIdx) - 0.00001
                % display feedback at the end of the block
                if target.catchTrials && ~isempty(strfind(targetBlockName,'abs'))
                    fixSeq(iFrame,1) = 9; % gray
                else
                    fixSeq(iFrame,1) = 8; % blue
                end
            else
                fixSeq(iFrame,1) = 5;
            end
        otherwise
            error('attBlockName not recognized')
    end
    
    % determine correct response keycode
    cueBlock = cueBlockNames{cueBlockOrder(blockIdx)};
    if ~strcmp(cueBlock, 'no-cue') && ...
            time-blockStartTimes(blockIdx) > targetLeadTime + jitSeq(blockIdx) + targetSOA + cueTargetSOA && ...
            blockStartTimes(blockIdx+1)-time > feedbackDur + itiSeq(blockIdx) - 0.00001
        % which target is post-cued?
        responseCue = str2double(cueBlock(end));
        switch targetBlockName
            case 'pres-pres'
                targetFrames = targetTypeSeq(find(targetTypeSeq>0,...
                    nFramesPerTarget*2,'last'));
                targets = targetFrames([1 end])';
                
                posFrames = posSeq(find(targetTypeSeq>0,...
                    nFramesPerTarget*2,'last'));
                positions = posFrames([1 end])';
                
                pedFrames = pedestalSeq(find(targetTypeSeq>0,...
                    nFramesPerTarget*2,'last'));
                pedestals = pedFrames([1 end])';
                
                switch responseOption
                    case 'targetType'
                        correctResponse = targets(responseCue);
                    case 'targetPos'
                        [row, col] = find(target.responsePosSets==positions(responseCue));
                        correctResponse = col;
                    case 'targetContrast4Levels'
                        correctResponse = 2*(pedestals(responseCue)-1) + targets(responseCue);
                    otherwise
                        error('responseOption not recognized')
                end
            case 'pres-abs'
                targetFrames = targetTypeSeq(find(targetTypeSeq>0,...
                    nFramesPerTarget,'last'));
                targets = [targetFrames(1) 0];
                
                posFrames = posSeq(find(targetTypeSeq>0,...
                    nFramesPerTarget,'last'));
                positions = [posFrames(1) 0];
                
                pedFrames = pedestalSeq(find(targetTypeSeq>0,...
                    nFramesPerTarget,'last'));
                pedestals = [pedFrames(1) 0];
                
                if responseCue==1
                    switch responseOption
                        case 'targetType'
                            correctResponse = targets(1);
                        case 'targetPos'
                            [row col] = find(target.responsePosSets==positions(1));
                            correctResponse = col;
                        case 'targetContrast4Levels'
                            correctResponse = 2*(pedestals(1)-1) + targets(1);
                    end
                else
%                     if target.catchTrials
%                         % random correct response for absent catch trials
%                         correctResponse = round(rand)+1;
%                     else
                        correctResponse = 3; % absent
%                     end
                end
            case 'abs-pres'
                targetFrames = targetTypeSeq(find(targetTypeSeq>0,...
                    nFramesPerTarget,'last'));
                targets = [0 targetFrames(1)];
                
                posFrames = posSeq(find(targetTypeSeq>0,...
                    nFramesPerTarget,'last'));
                positions = [0 posFrames(1)];
                
                pedFrames = pedestalSeq(find(targetTypeSeq>0,...
                    nFramesPerTarget,'last'));
                pedestals = [0 pedFrames(1)];
                
                if responseCue==1
%                     if target.catchTrials
%                         % random correct response for absent catch trials
%                         correctResponse = round(rand)+1;
%                     else
                        correctResponse = 3; % absent
%                     end
                else
                    switch responseOption
                        case 'targetType'
                            correctResponse = targets(2);
                        case 'targetPos'
                            [row, col] = find(target.responsePosSets==positions(2));
                            correctResponse = col;
                        case 'targetContrast4Levels'
                            correctResponse = 2*(pedestals(2)-1) + targets(2);
                    end
                end
            case 'abs-abs'
                targets = [0 0];
                positions = [0 0];
                correctResponse = 3; % absent
        end
        trialsPresented(blockIdx).precue = str2double(cueBlockName(1));
        trialsPresented(blockIdx).responseCue = responseCue;
        trialsPresented(blockIdx).targets = targets;
        trialsPresented(blockIdx).positions = positions;
        trialsPresented(blockIdx).pedestals = pedestals;
        trialsPresented(blockIdx).correctResponse = correctResponse;
        keyCodeSeq(iFrame,1) = keyCodes(correctResponse);
    else
        keyCodeSeq(iFrame,1) = 0;
    end
    
    switch triggerOption
        case 'conditionID'
            % determine trigger - condition ID
            if newBlock % only give condition trigger at the first frame of the block
                switch cueBlock
                    case 'no-cue' % blank
                        trig = 7;
                    case '1-1'
                        trig = 1;
                    case '1-2'
                        trig = 2;
                    case '2-1'
                        trig = 3;
                    case '2-2'
                        trig = 4;
                    otherwise
                        error('cueBlock not recognized')
                end
            elseif targetOnSeq(iFrame)~=0
                % triger for target side, only on first target frame
                if targetOnSeq(iFrame)==1 && targetOnSeq(iFrame-1)==0
                    trig = NaN; % target on left % never happens, so don't use up the trigger
                elseif targetOnSeq(iFrame)==2 && targetOnSeq(iFrame-1)==0
                    trig = 6; % target on right
                else
                    trig = NaN;
                end
            elseif targetAbsOnSeq(iFrame)~=0
                % triger for target side, only on first target frame
                if targetAbsOnSeq(iFrame)==1 && targetAbsOnSeq(iFrame-1)==0
                    trig = NaN; % target on left % never happens, so don't use up the trigger
                elseif targetAbsOnSeq(iFrame)==2 && targetAbsOnSeq(iFrame-1)==0
                    trig = NaN; % 5; target on right; used to be 5 - now never happens 
                else
                    trig = NaN;
                end
            elseif cueSeq(iFrame)~=0
%                 trig = 8; % pre- or post-cue
                if prePostCueSeq(iFrame) == 1
                    trig = 8; % pre-cue
                else
                    trig = 5; % post-cue
                end
            else
                trig = NaN;
            end
            
            trigSeq(iFrame,1) = computeTrigger(trig);
            
        otherwise
            error('triggerOption not recognized')
    end
end

% % show targetOnSeq and trigSeq
% f(1) = figure;
% subplot(3,1,1)
% plot(seqtiming,targetOnSeq)
% subplot(3,1,2)
% plot(seqtiming,trigSeq)
% subplot(3,1,3)
% plot(seqtiming,keyCodeSeq)
% 
% % display triggers by channel
% f(2) = displayTrigger(trigSeq, nBlocks);
% set(f(2),'Position',[0 0 1200 900]);

%% Create stimulus strucutre
% set remaining stimulus variables
cmap = repmat((0:255)',1,3);
srcRect = [0 0 size(images,2) size(images,1)];
destRect = CenterRectOnPoint(srcRect, cx, cy); %(cy+stimPos(2)*pixelsPerDegree)
diodeSeq = repmat(slowUnit-1, 1, ceil(length(seq)/length(slowUnit)))';
target.seq = targetTypeSeq;
target.posSeq = posSeq;
target.pedestalSeq = pedestalSeq;

% store in stimulus structure
stimulus.fixDiam = fixDiam;
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
stimulus.target = target;
stimulus.respDur = respDur;
stimulus.itiSeq = itiSeq; % storage only
stimulus.jitSeq = jitSeq;

% friendlier format of posBlockOrder
posBlockOrder = nan(size(targetTypeBlockOrder));
for iBlock = 1:nBlocks-1
    positions = trialsPresented(iBlock).positions;
    if ~isempty(positions)
        posBlockOrder(:,iBlock) = positions;
    end
end

% friendlier format of pedestalBlockOrder
pedestalBlockOrder = nan(size(targetTypeBlockOrder));
for iBlock = 1:nBlocks-1
    pedestals = trialsPresented(iBlock).pedestals;
    if ~isempty(pedestals)
        pedestalBlockOrder(:,iBlock) = pedestals;
    end
end

% store in order structure
order.blockOrder = blockOrder;
order.attBlockOrder = attBlockOrder;
order.targetBlockOrder = targetBlockOrder;
order.cueBlockOrder = cueBlockOrder; 
order.targetTypeBlockOrder = targetTypeBlockOrder;

if strcmp(target.type, 'contrast')
    order.posShuffled = posShuffled; %8 positions for each condition (pres-pres, pres-abs, abs-pres)
    order.posBlockOrder = posBlockOrder;
%     order.pedestalShuffled = pedestalShuffled;
    order.targetPedestals = targetPedestals;
    order.pedestalBlockOrder = pedestalBlockOrder;
end
order.targetTypes = targetTypes;
order.trialsPresented = trialsPresented;

% save stimulus
if saveStim
    save(sprintf('%s/%s.mat', stimDir, stimFile), 'stimulus', 'p','order')
    save(sprintf('%s/%s.mat', vistaStimDir, stimFile), 'stimulus', 'p','order')
end

% save figs
if saveFigs
    rd_saveAllFigs(f, {'trigplot','trigchan'}, stimFile);
end
