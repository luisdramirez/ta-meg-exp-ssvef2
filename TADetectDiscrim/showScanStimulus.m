function [response, timing, quitProg] = showScanStimulus(display,...
    stimulus, t0, timeFromT0)
% [response, timing, quitProg] = showStimulus(display,stimulus, ...
%           [time0 = GetSecs], [timeFromT0 = true])
%
% Inputs
%   display:    vistaisp display structure
%   stimulus:   vistadisp stimulus structure (e.g., see doRetinotopyScan.m)
%   t0:         time the scan started in seconds acc to PsychtoolBox
%               GetSecs function. By default stimulus timing is relative to
%               t0. If t0 does not exist it is created at the start of this
%               program.
%   timeFromT0: boolean. If true (default), then time each screen flip from
%               t0. If false, then time each screen flip from last screen
%               flip. The former is typically used for fMRI, where we want
%               to avoid accumulation of timing errors. The latter may be
%               more useful for ECoG/EEG where we care about the precise
%               temporal frequency of the stimulus.
% Outputs:
%   response:   struct containing fields
%                   keyCode: keyboard response at each frame, if any; if
%                           no response record a 0);
%                   secs: time of each response in seconds ?? verify
%                   flip:   time of each screen flip measured by PTB
%   timing:     float indicating total time of experiment
%   quitProg:   Boolean to indicate if experiment ended by hitting quit key
%
%
% HISTORY:
% 2005.02.23 RFD: ported from showStimulus.
% 2005.06.15 SOD: modified for OSX. Use internal clock for timing rather
%                 than framesyncing because getting framerate does not
%                 always work. Using the internal clock will also allow
%                 some "catching up" if stimulus is delayed for whatever
%                 reason. Loading mex functions is slow, so this should be
%                 done before callling this program.
% 2011.09.15  JW: added optional input flag, timeFromT0 (default = true).
%                 true, we time each screen flip from initial time (t0). If
%                 false, we time each screen flip from the last screen
%                 flip. Ideally the results are the same.

% flip every n frames
nFramesPerFlip = 3;

% triggers?
triggersOn = 1;

% staircase? (adjustment between runs)
staircase = 0; % set to 0 for first run of the day

% set target difficulty
if staircase && exist('staircase.mat','file')
    s = load('staircase.mat');
    switch stimulus.target.type
        case 'cb'
            tilts = [-s.tilt s.tilt];
            patchContrast = s.contrast;
            fprintf('\n\n\n\nSTAIRCASE UPDATE:\ncontrast = %.2f\ntilt = %1.1f\n\n\n\n', patchContrast, tilts(2))
        case 'contrast'
            patchContrast = s.contrasts;
            fprintf('\n\n\n\nSTAIRCASE UPDATE:\ncontrast = [%1.2f %1.2f]\n\n\n\n', patchContrast(1), patchContrast(2))
        otherwise
            error('stimulus.target.type not recognized')
    end
else
    % MANUAL SETTINGS
    tilts = [-2.5 2.5]; % starting settings: [-6 6] [relative to the base orientation]
    patchContrast = 0.4;
%     patchContrast = [.1 .2 .85 .95]; % rd-40
%     patchContrast = [0 .2 .8 .95]; % starting values 40
    dotSize = 0.3; % in degrees
    shifts = [0 0]; % phase shifts
    % patchSize = 1; % for cb target (this should be set in makeTADetectDiscrimStim, but just testing for now)
end
if strfind(display.position, 'MEG')
    soundAmp = 0.10; % 0.10 for MEG
else
    soundAmp = 1;
end

% input checks
if nargin < 2
    help(mfilename);
    return;
end;
if nargin < 3 || isempty(t0)
    t0 = GetSecs; % "time 0" to keep timing going
end;

if notDefined('timeFromT0'), timeFromT0 = true; end

% kbQueue
if ~isfield(display.devices,'useKbQueue')
    useKbQueue = 0;
else
    useKbQueue = display.devices.useKbQueue;
end

% some more checks
if ~isfield(stimulus,'textures')
    % Generate textures for each image
    disp('WARNING: Creating textures before stimulus presentation.');
    disp(['         This should be done before calling ' mfilename ' for']);
    disp('         accurate timing.  See "makeTextures" for help.');
    stimulus = makeTextures(display,stimulus);
end;

% quit key
if checkfields(display, 'quitProgKey')
    quitProgKey = display.quitProgKey;
else
    quitProgKey = KbName('q');
end;

% some variables
nFrames = length(stimulus.seq);
HideCursor;
nGamma = size(stimulus.cmap,3);
nImages = length(stimulus.textures);
response.keyCode = zeros(length(stimulus.seq),1); % get 1 buttons max
response.secs = zeros(size(stimulus.seq));        % timing
response.RT = zeros(size(stimulus.seq));
quitProg = 0;
response.flip = [];
lastSoundTime = NaN;

% set up KbQueue if desired
if useKbQueue
    KbQueueCreate(display.devices.keyInputExternal);
    KbQueueStart();
end

% set up sound if desired
if isfield(stimulus, 'soundSeq')
    % Perform basic initialization of the sound driver
    InitializePsychSound(1); % 1 for precise timing
    
    % Open audio device for low-latency output
    Fs = 44100;
    reqlatencyclass = 2; % Level 2 means: Take full control over the audio device, even if this causes other sound applications to fail or shutdown.
    pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, Fs, 1); % 1 = single-channel
    
    % Play example sounds
    Screen('Flip',display.windowPtr);
    % Wait for key press
    if useKbQueue
        KbQueueWait();
    else
        KbWait(-1);
    end
    WaitSecs(1)
    for iSound = 1:size(stimulus.sounds,1)
        playSound(pahandle, stimulus.sounds(iSound,:)*soundAmp);
        WaitSecs(1)
    end
end

% set up target if desired
if isfield(stimulus, 'target')
    target = stimulus.target;
    
    switch target.type
        case 'lines'
            fprintf('\n[showScanStimulus] lines tilt = [%1.1f %1.1f]\n\n', tilts(1), tilts(2))
            target.tilts = tilts; % store settings
        case 'dot'
            fprintf('\n[showScanStimulus] dot size = %1.1f degrees\n\n', dotSize)
            target.dotSize = dotSize; % store settings
            dotSizePx = dotSize*target.pixelsPerDegree;
        case 'grating'
            fprintf('\n[showScanStimulus] grating tilt = [%1.1f %1.1f], shift = [%1.2f %1.2f]\n\n', tilts(1), tilts(2), shifts(1)/pi, shifts(2)/pi)
            target.tilts = tilts; % store settings
            target.shifts = shifts;
            for iShift = 1:numel(shifts) % here, shift refers to either a phase shift or an orientation change, depending on the values of "shifts" and "tilts"
                for iP1 = 1:numel(target.phases)
                    for iP2 = 1:numel(target.phases)
                        p2 = target.phases(iP2) + shifts(iShift);
                        o2 = target.orientation + tilts(iShift);
                        
                        targ0 = buildColorGrating(target.pixelsPerDegree, [target.stimSize target.stimSize], ...
                            target.spatialFreq, o2, p2, target.contrast, ...
                            1, 'bw', 1, 1);
                        
                        targ1 = maskWithAnnulus(targ0, length(targ0), ...
                            0, target.blurRadius, target.backgroundColor);
                        
                        targ{iP1, iP2, iShift} = [target.stim{iP1,1} target.spacer targ1];
                        
                        target.textures(iP1, iP2, iShift) = Screen('MakeTexture', ...
                            display.windowPtr, targ{iP1,iP2,iShift}.*255);
                    end
                end
            end
        case 'cb'
            fprintf('\n[showScanStimulus] cb tilt = [%1.1f %1.1f], contrast = %1.2f\n\n', tilts(1), tilts(2), patchContrast)
            target.tilts = tilts; % store settings
            target.contrast = patchContrast;
            % target.size = patchSize;
            if strcmp(target.stimType,'noise')
                % add to background, no transparency
                tt = target.targetTypes;
                tp = target.targetPedestals;
                ntrials = length(stimulus.itiSeq);
                blankidx = 1:5:ntrials;
                targidx = setdiff(1:ntrials, blankidx);
                for iT = 1:numel(tt)
                    trialNum = ceil(iT/2);
                    tilt = (tp(iT)-1)*90 + tilts(tt(iT));
                    for iPhase = 1:2
                        phase = target.phases(iPhase);
                        targ0 = buildColorGrating(target.pixelsPerDegree, [target.imSize target.imSize], ...
                            target.spatialFreq, tilt, phase, target.contrast, 0, 'bw');
                        bgim = stimulus.images(:,:,targidx(trialNum)*2-1+iPhase-1);
                        if size(bgim,1)~=size(targ0,1)
                            if size(targ0,1)==size(bgim,1)+2
                                targ0 = targ0(2:end-1,2:end-1);
                            else 
                                targ0 = targ0(2:end,2:end);
                            end
                        end
                        targ1 = (targ0-.5) + (bgim/255-.5) + .5;
                        targ = maskWithAnnulus(targ1, size(bgim,1), 0, ...
                            target.blurRadius, target.backgroundColor);
                        target.textures(iT,iPhase) = Screen('MakeTexture', display.windowPtr, targ.*255);
                    end
                end
            else
                % use transparency
                Screen('BlendFunction', display.windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                
                targ0 = buildColorGrating(target.pixelsPerDegree, [target.imSize target.imSize], ...
                    target.spatialFreq, 0, 0, 1, 0, 'bw');
                alphaLayer = make2DGaussianCentered(size(targ0,1), size(targ0,2), 0, 0, target.size*target.pixelsPerDegree, 1).*patchContrast; % black is transparent, white is opaque
                targ = cat(3, targ0, alphaLayer);
                
                target.textures = Screen('MakeTexture', display.windowPtr, targ.*255);
            end
            targetSizePx = size(targ0);
            target.destRect(1:2) = target.center - targetSizePx/2;
            target.destRect(3:4) = target.destRect(1:2) + targetSizePx;
            target.baseOrients = []; % initialize
        case 'contrast'
            % CREATE TARGET IMAGES HERE based on parameters specified in
            % 'makeTADetectDiscrimStim'
%s            fprintf('\n[showScanStimulus] contrast: contrasts = [%1.2f %1.2f]\n\n', patchContrast(1), patchContrast(2))
            target.contrast = patchContrast; 
            [backgroundIms, maskedIms, target] = contrastTarget(display, target);
            [madeBgIms, madeTargIms] = makeIms(display, backgroundIms, maskedIms);
        otherwise
            error('target.type not recognized')
    end
end

% respose duration
if isfield(stimulus, 'respDur')
    respDur = stimulus.respDur;
    nRespFrames = display.frameRate*respDur; % look respDur s back
end
startedFeedback = 0;
targetCounter = 0;

% go
fprintf('[%s]:Running. Hit %s to quit.\n',mfilename,KbName(quitProgKey));

% If we are doing ECoG, then start with black photodiode
if isfield(stimulus, 'diodeSeq'), drawTrig(display,0); end

for frame = 1:nFramesPerFlip:nFrames
    
    %--- update display
    % If the sequence number is positive, draw the stimulus and the
    % fixation.  If the sequence number is negative, draw only the
    % fixation.
    if stimulus.seq(frame)>0
        % put in an image
        imgNum = mod(stimulus.seq(frame)-1,nImages)+1;
        Screen('DrawTexture', display.windowPtr, stimulus.textures(imgNum), stimulus.srcRect, stimulus.destRect);

        % add target
        if isfield(stimulus, 'target')
            if target.seq(frame)>0
                switch target.type
                    case 'lines'
                        rot = target.baseOrient + tilts(target.seq(frame));
                        xy = rotateCoords(target.xy0, rot);
                        Screen('DrawLines', display.windowPtr, xy, target.width, ...
                            target.colors, target.center);
                    case 'dot'
                        if target.seq(frame-1)==0 % update only for new targets
                            r = target.maxRadiusPx*rand;
                            theta = 2*pi*rand;
                            [x y] = pol2cart(theta, r);
                            y = abs(y)*target.dotLocs(target.seq(frame)); % pos or neg y
                            rect(1:2) = target.center + [x y] - dotSizePx/2;
                            rect(3:4) = rect(1:2) + dotSizePx;
                        end
                        Screen('FillOval', display.windowPtr, target.colors, rect);
                    case 'grating'
                        ph = target.phSeq(frame,:); % phase
                        sh = target.seq(frame); % shift
                        Screen('DrawTexture', display.windowPtr, target.textures(ph(1),ph(2),sh), ...
                            stimulus.srcRect, stimulus.destRect);
                    case 'cb'
                        if strcmp(target.stimType,'noise')
                             if target.seq(frame-1)==0 % update counter for new targets
                                targetCounter = targetCounter+1;
                                baseOrient = (target.targetPedestals(targetCounter)-1)*90;
                                target.baseOrients = [target.baseOrients baseOrient]; % store a copy
                             end
                             phaseIdx = 2-mod(stimulus.seq(frame),2);
                             Screen('DrawTexture', display.windowPtr, target.textures(targetCounter,phaseIdx), ...
                                 [], target.destRect);
                        else
                            if target.seq(frame-1)==0 % update only for new targets
                                baseOrient = round(rand)*90; % 0 or 90
                                target.baseOrients = [target.baseOrients baseOrient]; % store
                            end
                            rot  = baseOrient + tilts(target.seq(frame));
                            Screen('DrawTexture', display.windowPtr, target.textures(1), ...
                                [], target.destRect, rot);
                        end
                    case 'contrast'
                        % DRAW TARGET HERE   
                        if isfield(target,'pedestalSeq')
                            targetIdx = (target.pedestalSeq(frame)-1)*2 + target.seq(frame);
                        else
                            targetIdx = target.seq(frame);
                        end
                        if stimulus.seq(frame) == 3 
                            %bgtex = Screen('MakeTexture', display.windowPtr, backgroundIms{1}*255); OR madeBgIms{1}; 
                            %tex = Screen('MakeTexture', display.windowPtr, maskedIms{target.posSeq(frame), 1, target.seq(frame)}*255); OR madeTargIms{target.posSeq(frame),1,target.seq(frame)}; 
                            Screen('DrawTexture', display.windowPtr, madeBgIms{1});
                            Screen('DrawTexture', display.windowPtr, madeTargIms{target.posSeq(frame),1,targetIdx});
                        elseif stimulus.seq(frame) == 4 
                            %bgtex = Screen('MakeTexture', display.windowPtr, backgroundIms{2}*255); OR madeBgIms{2}; 
                            %tex = Screen('MakeTexture', display.windowPtr, maskedIms{target.posSeq(frame), 2, target.seq(frame)}*255); OR madeTargIms{target.posSeq(frame), 2, target.seq(frame)};          
                            Screen('DrawTexture', display.windowPtr, madeBgIms{2});
                            Screen('DrawTexture', display.windowPtr, madeTargIms{target.posSeq(frame), 2, targetIdx});
                        end
                    otherwise
                        error('target.type not recognized');
                end
            end
        end
        
        % determine feedback according to response accuracy
        % sorry some of this is hard-coded for now
        if stimulus.fixSeq(frame)>=8 % feedback period
            % if first frame of feedback period, determine accuracy
            if stimulus.fixSeq(frame-nFramesPerFlip)<8 && ~startedFeedback
                startedFeedback = 1;
                responseWindow = response.correct(frame-nRespFrames:frame-nFramesPerFlip);
                correct = responseWindow(responseWindow~=0);
            elseif frame==nFrames || stimulus.fixSeq(frame+1)<8 % last feedback frame
                startedFeedback = 0;
            end
            % change fixation according to accuracy
            if ~isempty(correct) && stimulus.fixSeq(frame)~=9 % keep 9 (catch trial) regardless of response
                correct = correct(1); % take first response
                if correct==1
                    stimulus.fixSeq(frame) = 7; % green
                elseif correct==-1
                    stimulus.fixSeq(frame) = 6; % red
                end
            end
        end
        
        %drawFixation_rd(display,stimulus.fixSeq(frame)); % FIXATION DRAWN HERE
        contrastFixation(display, stimulus, frame)
        
        % If we are doing eCOG, then flash photodiode if requested
        if isfield(stimulus, 'diodeSeq')
            colIndex = drawTrig(display,stimulus.diodeSeq(frame));
        end
        
        % Play a sound if requested
        if isfield(stimulus, 'soundSeq')
            if stimulus.soundSeq(frame)~=0
                playSound(pahandle, stimulus.sounds(stimulus.soundSeq(frame),:)*soundAmp);
                lastSoundTime = GetSecs;
            end
        end
        
    elseif stimulus.seq(frame)<0
        % put in a color table
        gammaNum = mod(-stimulus.seq(frame)-1,nGamma)+1;
        % The second argument is the color index.  This apparently changed
        % in recent times (07.14.2008). So, for now we set it to 1.  It may
        % be that this hsould be
        drawFixation_rd(display,stimulus.fixSeq(frame));
        Screen('LoadNormalizedGammaTable', display.windowPtr, stimulus.cmap(:,:,gammaNum));
    end;
    
    %--- timing
    [waitTime, requestedFlipTime] = getWaitTime(stimulus, response, frame,  t0, timeFromT0, nFramesPerFlip);
    
    %--- get inputs (subject or experimentor)
    % check if this is the last frame of the response window
    while(waitTime<0)
        % Scan the keyboard for subject response
        if useKbQueue
            % Use KbQueue
            [keyIsDown, firstPress] = KbQueueCheck();
            if keyIsDown
                secs = min(firstPress(firstPress~=0));
                ssKeyCode = firstPress==secs;
                response.keyCode(frame) = find(ssKeyCode);
                response.secs(frame) = secs - t0;
                response.RT(frame) = secs - lastSoundTime;
                
                if response.keyCode(frame)==stimulus.keyCodeSeq(frame)
                    response.correct(frame) = 1;
                else
                    response.correct(frame) = -1;
                end
                
                if(ssKeyCode(quitProgKey))
                    quitProg = 1;
                    break; % out of while loop
                end;
                if(response.correct(frame)~=0) %%% TO ADD OTHERWISE OVERWRITE THE response.correct VARIABLE
                    break; % out of while loop
                end;
            else
                response.correct(frame) = 0;
                if(response.correct(frame)==0) %%% TO ADD OTHERWISE OVERWRITE THE response.correct VARIABLE
                    break; % out of while loop
                end;
            end
        else
            % Use KbCheck
            %[ssKeyIsDown,ssSecs,ssKeyCode] = KbCheck(display.devices.keyInputExternal);
            [ssKeyIsDown,ssSecs,ssKeyCode] = KbCheck(-1);
            if(ssKeyIsDown)
                kc = find(ssKeyCode);
                response.keyCode(frame) = kc(1);
                % response.keyCode(frame) = 1; % binary response for now
                response.secs(frame)    = ssSecs - t0;
                response.RT(frame) = ssSecs - lastSoundTime;
                
                if response.keyCode(frame)==stimulus.keyCodeSeq(frame)
                    response.correct(frame) = 1;
                else
                    response.correct(frame) = -1;
                end
                
                if(ssKeyCode(quitProgKey))
                    quitProg = 1;
                    break; % out of while loop
                end;
            else
                response.correct(frame) = 0;
            end
        end
        
        % if there is time release cpu
        if(waitTime<-0.03), WaitSecs(0.01); end;
        
        
        % timing
        [waitTime, requestedFlipTime] = getWaitTime(stimulus, response, frame, t0, timeFromT0, nFramesPerFlip);
    end;
    
    %--- stop?
    if quitProg
        fprintf('[%s]:Quit signal recieved.\n',mfilename);
        break;
    end;
    
    %--- update screen
    VBLTimestamp = Screen('Flip', display.windowPtr, requestedFlipTime);

    % send trigger for MEG, if requested, and record the color of the PD
    % cue
    if isfield(stimulus, 'trigSeq') && stimulus.trigSeq(frame) > 0 && triggersOn
        PTBSendTrigger(stimulus.trigSeq(frame), 0);
%         fprintf('Trigger %d\n', stimulus.trigSeq(frame)); drawnow
        response.trig(frame) = stimulus.trigSeq(frame);
    end
     
    if isfield(stimulus, 'diodeSeq')
        response.LED(frame)  = colIndex;
    end
    
    % record the flip time
    response.flip(frame) = VBLTimestamp;
    
end;

% store target data inside response, so it gets passed out
target.soundAmp = soundAmp;
response.target = target;

% clean up KbQueue
if useKbQueue
    KbQueueStop();
    KbQueueRelease();
end

% that's it
ShowCursor;
timing = GetSecs-t0;
fprintf('[%s]:Stimulus run time: %f seconds [should be: %f].\n',mfilename,timing,max(stimulus.seqtiming));

return;


function [waitTime, requestedFlipTime] = getWaitTime(stimulus, response, frame, t0, timeFromT0, nFramesPerFlip)
% [waitTime, requestedFlipTime] = getWaitTime(stimulus, response, frame, t0, timeFromT0, nFramesPerFlip)
%
% If timeFromT0 we wait until the current time minus the initial time is
% equal to the desired presentation time, and then flip the screen.
% If timeFromT0 is false, then we wait until the current time minus the
% last screen flip time is equal to the desired difference in the
% presentation time of the current flip and the prior flip.

if timeFromT0
    waitTime = (GetSecs-t0)-stimulus.seqtiming(frame);
    requestedFlipTime = t0 + stimulus.seqtiming(frame);
else
    if frame > nFramesPerFlip
        lastFlip = response.flip(frame-nFramesPerFlip);
        desiredWaitTime = stimulus.seqtiming(frame) - stimulus.seqtiming(frame-nFramesPerFlip);
    else
        lastFlip = t0;
        desiredWaitTime = stimulus.seqtiming(frame);
    end
    % we add 10 ms of slop time, otherwise we might be a frame late.
    % This should NOT cause us to be 10 ms early, because PTB waits
    % until the next screen flip. However, if the refresh rate of the
    % monitor is greater than 100 Hz, this might make you a frame
    % early. [So consider going to down to 5 ms? What is the minimum we
    % need to ensure that we are not a frame late?]
    waitTime = (GetSecs-lastFlip)-desiredWaitTime + .030;
    requestedFlipTime = lastFlip + desiredWaitTime - 0.008;
end
