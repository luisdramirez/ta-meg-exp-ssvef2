function [responseData, runLabels] = sl_responseData(respTime, ...
    trialCount, respSecs,feedbackDur,refreshRate, blockLength, keyCodes, response, ...
    order,runNum)

% use like this: 
% load('taDetect7','order');        % load block order from stimFile
% load('20141028T144433_taDetect5.mat') % load data
% trialCount = 41;  
% respSecs = 1.4;
% feedbackDur = 0.3;
% refreshRate = 60;  %(frames)
% blockLength = 300; %(frames)
% respTime = 198;  % frames to respond period
% keyCodes = [30 31];
% runNum = 1;


for trialNum = 1 : trialCount
    
    % Defining search range
    begSearchIdx = (blockLength * (trialNum - 1) + respTime) + 1;
    endSearchIdx = begSearchIdx + ((respSecs + feedbackDur) * refreshRate - 1);
    
    % Find the index for ALL non-zero elements (key-codes) within search
    % range
    
    foundIdxArr = begSearchIdx + find(response.keyCode(begSearchIdx:endSearchIdx)) - 1;
    
    % ... If nothing is found let me know
    if length(foundIdxArr) < 1
        sprintf('Trial %d: No KeyCode value found', trialNum);
        % If no response, make these variables blank
        rData.keyCode(trialNum,:) = NaN;
        rData.RT(trialNum,:) = NaN;
        rData.response(trialNum,:) = NaN;
        rData.correct(trialNum,:) = NaN;
        % Since nothing was found just skip over to the next trial
        continue
    end
    
    % This is the index of the first key-code found
    fidx = foundIdxArr(1);
    % Use this to find the key code value
    keyVal = response.keyCode(fidx);
    RespVal = response.correct(fidx);
    % Prints info about the key code value found
    sprintf('Trial %d: KeyCode %d found at Idx = %d', trialNum, keyVal, fidx);
   
    % Stores key code found in data matrix
    rData.RT(trialNum,:) = response.secs(fidx) - ((blockLength * (trialNum - 1) + respTime)/refreshRate);
    rData.keyCode(trialNum,:) = keyVal; 
    rData.correct(trialNum,:) = RespVal;
    
    % Create translation column for what keyCode means (response)
    switch keyVal
        case keyCodes(1)
            rData.response(trialNum,:) = 1;
        case keyCodes(2)
            rData.response(trialNum,:) = 2;
    end
end



%% Create final output matrix

runLabels = {'run number', 'trial number','fast side','cue condition','attended side',... 
    'target condition','key code', 'response', 'correct', 'RT'};
responseData = horzcat(repmat(runNum,trialCount,1),(1:trialCount)', order.blockorder', ...
    order.cueBlockOrder',order.attBlockOrder',order.targetBlockOrder', ...
    rData.keyCode, rData.response,rData.correct, rData.RT);

end


