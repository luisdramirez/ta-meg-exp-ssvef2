function [blockOrder, attBlockOrder, targetBlockOrder, cueBlockOrder, targetTypeBlockOrder, targetPedestalBlockOrder] = ...
    block_gen3(subjectID, run, stimDir)

% function [blockOrder, attBlockOrder, targetBlockOrder, cueBlockOrder, targetTypeBlockOrder, targetPedestalBlockOrder] = ...
%    block_gen3(subjectID, run, stimDir)
%
% Assumes that runs will be executed IN ORDER, i.e. 1, 2, 3, 4.
% If a run is not found, the code will automatically generate a set of 4
% consecutive runs starting with the input run.

if nargin==0
    subjectID = 'test';
    run = 1;
    stimDir = 'stimuli';
end

%% file i/o
saveDir = sprintf('%s/%s', stimDir, subjectID);
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

%% check to see if runs have already been generated
fileName = sprintf('%s/%s_blockgen_run%d.mat', saveDir, subjectID, run);
if exist(fileName,'file')
    generateRuns = 0;
else
    generateRuns = 1;
end

%% generate runs if needed
if generateRuns
    %% setup
    nContrastT1 = 4;
    nContrastT2 = 4;
    nV = 4;
    nT = 2;
    
    nTrialsPerMiniblock = 4;
    nRunsForCounterbalance = 4;
    nTrialsPerRun = 41;
    
    runs = run:run+nRunsForCounterbalance-1;
    
    %% counterbalance factors
    a = fullfact([nContrastT1 nContrastT2 nV nT]);
    
    validity = a(:,3);
    validity(validity~=4) = 1;
    validity(validity==4) = 2;
    a(:,3) = validity;
    
    %% randomize trial sequence
    idx = randperm(size(a,1));
%     idx = 1:size(a,1);
    b = a(idx,:);
    
    %% add blank trials
    blankIdx0 = [repmat(1:nTrialsPerMiniblock+1, 1, size(a,1)/(nRunsForCounterbalance*nTrialsPerMiniblock)) 1];
    blankIdx0(blankIdx0~=1) = 0;
    blankIdx = repmat(blankIdx0, 1, nRunsForCounterbalance);
    
    c = zeros(length(blankIdx), size(b,2));
    c(blankIdx==0,:) = b;
    
    nullOrder = blankIdx;
    nullOrder(blankIdx==0) = 2;
    
    %% read out conditions
    contrastT1 = c(:,1);
    contrastT2 = c(:,2);
    validity = c(:,3);
    target = c(:,4);
    
    [t1TT, t1TP] = contrastToTypePedestal(contrastT1);
    [t2TT, t2TP] = contrastToTypePedestal(contrastT2);
    
    targetType = [t1TT' t2TT'];
    targetPedestal = [t1TP' t2TP'];
    
    cue = zeros(size(target));
    cue(validity==1) = target(validity==1);
    cue(validity==2) = 3-target(validity==2);
    
    cueType = zeros(size(cue));
    cueType(cue==0) = 1;
    cueType(cue==1 & target==1) = 2;
    cueType(cue==1 & target==2) = 3;
    cueType(cue==2 & target==1) = 4;
    cueType(cue==2 & target==2) = 5;
    
    %% split into runs and save needed variables
    for iRun = 1:numel(runs)
        idx = (1:nTrialsPerRun) + nTrialsPerRun*(iRun-1);
        
        blockOrder = nullOrder(idx);
        attBlockOrder = nullOrder(idx);
        targetBlockOrder = nullOrder(idx);
        
        cueBlockOrder = cueType(idx)';
        targetTypeBlockOrder = targetType(idx,:)';
        targetPedestalBlockOrder = targetPedestal(idx,:)';
        
        saveFileName = sprintf('%s/%s_blockgen_run%d.mat', saveDir, subjectID, runs(iRun));
        save(saveFileName, 'blockOrder', 'attBlockOrder', 'targetBlockOrder', ...
            'cueBlockOrder', 'targetTypeBlockOrder', 'targetPedestalBlockOrder')
    end
end

%% load run to return variables as output
% clear variables first to be sure we output the correct ones
[blockOrder, attBlockOrder, targetBlockOrder, cueBlockOrder, targetTypeBlockOrder, targetPedestalBlockOrder] = deal([]);

load(fileName)

