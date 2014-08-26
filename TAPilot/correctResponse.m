% Determine correct/incorrect responses for MEG data.

%% To do
% Replace hard-coded numbers with named variables
% Make trial headers cell array
% What happens if someone misses a trial (doesn't press any key)?

%% Just get the data from one file
load('20140806T103117_taPilot7.mat');
trialCount = 8;
%keyIndex = ([1:8]*840)+1
% attendIndex = ([1:8]*420)+1
respSecs = 2;
refreshRate = 60;

%% Key-code response search loop
for trialNum = 1 : trialCount
    % Search from begSearchIdx to endSearchIdx
    begSearchIdx = (840 * trialNum) + 1;
    endSearchIdx = begSearchIdx + (respSecs * refreshRate);
    
    % Find the index for ALL non-zero key-codes
    foundIdxArr = begSearchIdx + find(response.keyCode(begSearchIdx:endSearchIdx)) - 1;
    % Did we find anything?
    if length(foundIdxArr) < 1
        % Didn't find activity ... just skip over to the next trial
        sprintf('Trial %d: No KeyCode value found', trialNum)
        continue
    end
    
    % This is the first key-code found
    fidx = foundIdxArr(1);
    keyVal = response.keyCode(fidx);
    % We DID find a KeyCode and it's here ...
    fprintf('\nTrial %d: KeyCode %d found at Idx = %d', trialNum, keyVal, fidx))
    
    % Your key-code related code goes here
    rData.RT(trialNum,:) = response.secs(fidx) - (840 * trialNum +1)/60 ;
    rData.keyCode(trialNum,:) = keyVal;
    switch keyVal
        case 30
            rData.response(trialNum,:) = 0;
        case 31
            rData.response(trialNum,:) = 1;
        case 32
            rData.response(trialNum,:) = 2;
        case 33
            rData.response(trialNum,:) = 3;
    end
end

%% The cue search loop
% for trialNum = 1 : trialCount
%     cueIdx = 421 + (trialNum - 1) * 840;
cond = response.trig~=0 & response.trig<10;
cueVal = response.trig(cond);
cueIdx = find(cond);
% We DID find a cue and it's here ...
fprintf('\nTrial %d: cue = %d found at Idx = %d', trialNum, cueVal, cueIdx)

% Your key-code related code goes here
rData.cueVal = cueVal';

for trialNum = 1:length(rData.cueVal)
    switch rData.cueVal(trialNum)
        case {1, 4}
            rData.cueSide(trialNum,:) = 1; % 1= left
        case {2, 8}
            rData.cueSide(trialNum,:) = 2; % 2 = right
    end
end
% end

%% Trigger type search loop
for trialNum = 1 : trialCount
    % Search from begSearchIdx to endSearchIdx
    begSearchIdx = 1 + 421 + (trialNum - 1) * 840;
    endSearchIdx = 840 * trialNum;
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
    sixteenCount = 0;
    thirtytwoCount = 0;
    for trigIdx = foundIdxArr
        trigType = response.trig(trigIdx);
        if trigType == 16
            sixteenCount = sixteenCount + 1;
        elseif trigType == 32
            thirtytwoCount = thirtytwoCount + 1;
        else
            sprintf('WARNING: What?  Bad trigger type value: %d', trigType)
        end
    end
    % Put them in an array for each trial
    switch rData.cueVal(trialNum)
        case {1, 4}
            rData.attendedStim(trialNum,:) = sixteenCount;
            rData.unattendedStim(trialNum,:) = thirtytwoCount;
        case {2, 8}
            rData.attendedStim(trialNum,:) = thirtytwoCount;
            rData.unattendedStim(trialNum,:) = sixteenCount;
    end
end

%% The accuracy loop
for trialNum = 1 : trialCount;
    if rData.response(trialNum) == rData.attendedStim(trialNum);
        rData.accuracy(trialNum,:) = 1;
    else
        rData.accuracy(trialNum,:) = 0;
    end
end

trial = horzcat( (1:8)' , rData.cueVal , rData.cueSide , rData.attendedStim , rData.unattendedStim , rData.keyCode , rData.response , rData.accuracy , rData.RT );

