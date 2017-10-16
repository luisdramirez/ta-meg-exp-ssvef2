% codetest.m

subjectID = 'test';

stimDir = sprintf('stimuli/%s', subjectID);

blockNames = {'blank','fast-left'}; % fast-left, slow-left
attBlockNames = {'no-att','att-right'}; % att-right
targetBlockNames = {'no-targ','pres-pres'};
cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'}; % 2-1 = cueT2,postcueT1

s = [];
for run = 1:4
    fileName = sprintf('%s/%s_blockgen_run%d.mat', stimDir, subjectID, run);
    load(fileName);
    s(run).blockOrder = blockOrder;
    s(run).attBlockOrder = attBlockOrder;
    s(run).targetBlockOrder = targetBlockOrder;
    s(run).cueBlockOrder = cueBlockOrder;
    s(run).targetTypeBlockOrder = targetTypeBlockOrder;
    s(run).targetPedestalBlockOrder = targetPedestalBlockOrder;
    
%     [s(run).blockOrder, s(run).attBlockOrder, s(run).targetBlockOrder, s(run).cueBlockOrder, s(run).targetTypeBlockOrder, s(run).targetPedestalBlockOrder] ...
%         = block_gen3(subjectID, run, stimDir);
end


trials = [];
nTrials = 41;
for iRun = 1:4
    for iTrial = 1:nTrials
        cueName = cueBlockNames{s(iRun).cueBlockOrder(iTrial)};
        switch cueName(1)
            case '1'
                cueTarget = 1;
            case '2'
                cueTarget = 2;
            otherwise
                cueTarget = NaN;
        end
        switch cueName(end)
            case '1'
                responseTarget = 1;
                responeTargetIdx = 1;
            case '2'
                responseTarget = 2;
                responseTargetIdx = 2;
            otherwise
                responseTarget = NaN;
                responseTargetIdx = 1;
        end
        trial = (iRun-1)*nTrials + iTrial;
        trials(trial,:) = [s(iRun).blockOrder(iTrial), ...
            s(iRun).attBlockOrder(iTrial), s(iRun).targetBlockOrder(iTrial), ...
            s(iRun).cueBlockOrder(iTrial), cueTarget, responseTarget, ...
            s(iRun).targetTypeBlockOrder(responseTargetIdx,iTrial), ...
            s(iRun).targetPedestalBlockOrder(responseTargetIdx,iTrial)];
    end
end

cueType = trials(:,4);
cueTarget = trials(:,5);
responseTarget = trials(:,6);
targetType = trials(:,7);
targetPedestal = trials(:,8);
targetContrast = (targetPedestal-1)*2 + targetType;
cueValidity = double(cueTarget==responseTarget);
cueValidity(isnan(cueTarget)) = NaN;

t = [cueValidity responseTarget targetContrast];

%%%%%
% cueValidity = a(:,3);
% cueValidity(cueValidity==2) = 0;
% responseTarget = a(:,4);
% targetContrast(a(:,4)==1,1) = a(a(:,4)==1,1);
% targetContrast(a(:,4)==2,1) = a(a(:,4)==2,2);
%%%%%

count = [];
vs = [0 1];
ts = [1 2];
cs = 1:4;
for iT = 1:2
    target = ts(iT);
    wT = responseTarget==target;
    for iV = 1:2
        validity = vs(iV);
        wV = cueValidity==validity;
        for iC = 1:4
            contrast = cs(iC);
            wC = targetContrast==contrast;
            count(iT,iV,iC) = nnz(wT & wV & wC);
        end
    end
end

