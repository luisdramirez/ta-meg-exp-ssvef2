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

% set target difficulty
tilts = [-5 5]; % relative to the base orientation
dotSize = 0.5; % in degrees

% input checks
if nargin < 2,
    help(mfilename);
    return;
end;
if nargin < 3 || isempty(t0),
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
quitProg = 0;
response.flip = [];

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
end

% set up target if desired
if isfield(stimulus, 'target')
    target = stimulus.target;
    dotSizePx = dotSize*target.pixelsPerDegree;
    
    switch target.type
        case 'lines'
            fprintf('\n[showScanStimulus] tilt = [%1.1f %1.1f]\n\n', tilts(1), tilts(2))
        case 'dot'
            fprintf('\n[showScanStimulus] dot size = %1.1f degrees\n\n', dotSize)
        otherwise
            error('target.type not recognized')
    end
end

% go
fprintf('[%s]:Running. Hit %s to quit.\n',mfilename,KbName(quitProgKey));

% If we are doing ECoG, then start with black photodiode
if isfield(stimulus, 'diodeSeq'), drawTrig(display,0); end

for frame = 1:nFrames
    
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
                    otherwise
                        error('target.type not recognized');
                end
            end
        end
        
        % determine feedback according to response accuracy
        % sorry some of this is hard-coded for now
        if stimulus.fixSeq(frame)==8 % feedback period
            nRespFrames = display.frameRate*1.7; % look 1.7 s back
            responseWindow = response.correct(frame-nRespFrames:frame-1);
            correct = responseWindow(responseWindow~=0);
            if ~isempty(correct)
                correct = correct(1); % take first response
                if correct==1
                    stimulus.fixSeq(frame) = 7; % green
                elseif correct==-1
                    stimulus.fixSeq(frame) = 6; % red
                end
            end
        end
        
        drawFixation_rd(display,stimulus.fixSeq(frame));
        
        % If we are doing eCOG, then flash photodiode if requested
        if isfield(stimulus, 'diodeSeq')
            colIndex = drawTrig(display,stimulus.diodeSeq(frame));
        end
        
        % Play a sound if requested
        if isfield(stimulus, 'soundSeq')
            if stimulus.soundSeq(frame)~=0
                playSound(pahandle, stimulus.sounds(stimulus.soundSeq(frame),:))
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
    waitTime = getWaitTime(stimulus, response, frame,  t0, timeFromT0);
    
    %--- get inputs (subject or experimentor)
    while(waitTime<0),
        % Scan the keyboard for subject response
        if useKbQueue
            % Use KbQueue
            [keyIsDown firstPress] = KbQueueCheck();
            if keyIsDown
                secs = min(firstPress(firstPress~=0));
                ssKeyCode = firstPress==secs;
                response.keyCode(frame) = find(ssKeyCode);
                response.secs(frame) = secs - t0;
                
                if response.keyCode(frame)==stimulus.keyCodeSeq(frame)
                    response.correct(frame) = 1;
                else
                    response.correct(frame) = -1;
                end
                
                if(ssKeyCode(quitProgKey)),
                    quitProg = 1;
                    break; % out of while loop
                end;
            else
                response.correct(frame) = 0;
            end
        else
            % Use KbCheck
            %[ssKeyIsDown,ssSecs,ssKeyCode] = KbCheck(display.devices.keyInputExternal);
            [ssKeyIsDown,ssSecs,ssKeyCode] = KbCheck;
            if(ssKeyIsDown)
                kc = find(ssKeyCode);
                response.keyCode(frame) = kc(1);
                % response.keyCode(frame) = 1; % binary response for now
                response.secs(frame)    = ssSecs - t0;
                
                if response.keyCode(frame)==stimulus.keyCodeSeq(frame)
                    response.correct(frame) = 1;
                else
                    response.correct(frame) = -1;
                end
                
                if(ssKeyCode(quitProgKey)),
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
        waitTime = getWaitTime(stimulus, response, frame, t0, timeFromT0);
        
    end;
    
    %--- stop?
    if quitProg,
        fprintf('[%s]:Quit signal recieved.\n',mfilename);
        break;
    end;
    
    %--- update screen
    VBLTimestamp = Screen('Flip',display.windowPtr);
    
    % send trigger for MEG, if requested, and record the color of the PD
    % cue
    if isfield(stimulus, 'trigSeq') && stimulus.trigSeq(frame) > 0
        PTBSendTrigger(stimulus.trigSeq(frame), 0);
        fprintf('Trigger sent, %s\n, %s', datestr(now), stimulus.trigSeq(frame)); drawnow
        response.trig(frame) = stimulus.trigSeq(frame);
    end
    
    if isfield(stimulus, 'diodeSeq') 
        response.LED(frame)  = colIndex;
    end
    
    % record the flip time
    response.flip(frame) = VBLTimestamp;
    
end;

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


function waitTime = getWaitTime(stimulus, response, frame, t0, timeFromT0)
% waitTime = getWaitTime(stimulus, response, frame, t0, timeFromT0)
%
% If timeFromT0 we wait until the current time minus the initial time is
% equal to the desired presentation time, and then flip the screen.
% If timeFromT0 is false, then we wait until the current time minus the
% last screen flip time is equal to the desired difference in the
% presentation time of the current flip and the prior flip.

if timeFromT0
    waitTime = (GetSecs-t0)-stimulus.seqtiming(frame);
else
    if frame > 1,
        lastFlip = response.flip(frame-1);
        desiredWaitTime = stimulus.seqtiming(frame) - stimulus.seqtiming(frame-1);
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
    waitTime = (GetSecs-lastFlip)-desiredWaitTime + .010;
end
