%% Determine correct/incorrect responses for MEG data.
    %% Revisions: how to deal with no reponse

%% Just get the data from one file
load('20140806T101303_taPilot1.mat');
trialCount = 8;
respSecs = 2;
refreshRate = 60;
blockLength = 840; % in frames
cueLeft = {1 4};
cueRight = {2 8};
maxCue = max(max(cell2mat(cueLeft)), max(cell2mat(cueRight)));
keyCode = [30, 31, 32, 33]; % for response button box
targetLeft = 16;
targetRight = 32;

%% Key-code response search loop

for trialNum = 1 : trialCount
    % Search from begSearchIdx to endSearchIdx
    begSearchIdx = (blockLength * trialNum) + 1;
    endSearchIdx = begSearchIdx + (respSecs * refreshRate);
    
    % Find the index for ALL non-zero key-codes
    foundIdxArr = begSearchIdx + find(response.keyCode(begSearchIdx:endSearchIdx)) - 1;
    % Did we find anything?
    if length(foundIdxArr) < 1
        % Didn't find activity ... just skip over to the next trial
        sprintf('Trial %d: No KeyCode value found', trialNum)
        rData.keyCode(trialNum,:) = NaN;
        rData.RT(trialNum,:) = NaN;
        rData.response(trialNum,:) = NaN;
        continue
    end

    % This is the first key-code found
    fidx = foundIdxArr(1);
    keyVal = response.keyCode(fidx);
    % We DID find a KeyCode and it's here ...
    sprintf('Trial %d: KeyCode %d found at Idx = %d', trialNum, keyVal, fidx)

    % Your key-code related code goes here
    rData.RT(trialNum,:) = response.secs(fidx) - (blockLength * trialNum +1)/refreshRate ;
    rData.keyCode(trialNum,:) = keyVal;
    
    % Create translation column for what keyCode means (response)
    switch keyVal
        case keyCode(1)
            rData.response(trialNum,:) = 0;
        case keyCode(2)
            rData.response(trialNum,:) = 1;
        case keyCode(3)
            rData.response(trialNum,:) = 2;
        case keyCode(4)
            rData.response(trialNum,:) = 3;
    end
end

%% The attentional cue search loop

    cond = response.trig~=0 & response.trig<=maxCue ; 
    cueVal = response.trig(cond);
    cueIdx = find(cond);
    % We DID find a cue and it's here ...
    fprintf('Trial %d: cue = %d found at Idx = %d', trialNum, cueVal, cueIdx);
    % Your key-code related code goes here
    rData.cueVal = cueVal';
    
    %Create translation column for what cueVal mean (cue side)
    for trialNum = 1:length(rData.cueVal)
        switch rData.cueVal(trialNum)
            case cueLeft
                rData.cueSide(trialNum,:) = 1; % 1= left
            case cueRight
                rData.cueSide(trialNum,:) = 2; % 2 = right
        end
    end


%% Trigger type search loop
for trialNum = 1 : trialCount
    
    % Search from begSearchIdx to endSearchIdx
    begSearchIdx = cueIdx(trialNum) +1; % +1 to start one further than cue trigger
    endSearchIdx = blockLength * trialNum;  % this is one before blank trigger
    
    % AVOID: Exceed matrix dimensions
    if endSearchIdx > length(response.trig)
        endSearchIdx = length(response.trig);
    end
    
    % NOTE: begSearchIdx and endSearchIdx would be adjusted specifically
    %       for response.trig
    foundIdxArr = begSearchIdx + find(response.trig(begSearchIdx:endSearchIdx)) - 1;
    % Did we find anything?
    if length(foundIdxArr) < 1
        % Didn't find activity ... just skip over to the next trial
        sprintf('Trial %d: No trigger values found', trialNum)
        continue
    end
    % Count the trigger types in different arrays
    targetLeftCount = 0;
    targetRightCount = 0;
    for trigIdx = foundIdxArr
        trigType = response.trig(trigIdx);
        if trigType == targetLeft
            targetLeftCount = targetLeftCount + 1;
        elseif trigType == targetRight
            targetRightCount = targetRightCount + 1;
        else
            sprintf('WARNING: What?  Bad trigger type value: %d', trigType)
        end
    end
    % Put them in an array for each trial
    switch rData.cueVal(trialNum)
        case cueLeft
            rData.attendedStim(trialNum,:) = targetLeftCount;
            rData.unattendedStim(trialNum,:) = targetRightCount;
        case cueRight
            rData.attendedStim(trialNum,:) = targetRightCount;
            rData.unattendedStim(trialNum,:) = targetLeftCount;
    end
end

%% The accuracy loop
for trialNum = 1 : trialCount;
    if rData.response(trialNum) == rData.attendedStim(trialNum);
        rData.accuracy(trialNum,:) = 1; %correct
    else
        rData.accuracy(trialNum,:) = 0; % incorrect
    end
end

%% Create final output matrix
trial = horzcat( (1:trialCount)' , rData.cueVal , rData.cueSide , rData.attendedStim , rData.unattendedStim , rData.keyCode , rData.response , rData.accuracy , rData.RT );

trialLabels = {'trial number', 'cue trigger', 'cued side', 'attended targets', 'unattended targets', 'key code', 'response', 'accuracy', 'RT'};
