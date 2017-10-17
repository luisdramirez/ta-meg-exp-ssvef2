function behav = behavior(behav)

%% key
% blockNames = {'blank','fast-left'}; % fast-left
% attBlockNames = {'no-att','att-right'}; % att-right
% targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres','abs-abs'};
% cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'}; % 2-1 = cueT2,postcueT1

%% organize data
nTrials = size(behav.responseData_all,1);

cueCondIdx = strcmp(behav.responseData_labels,'cue condition');
t1CondIdx = strcmp(behav.responseData_labels,'target type T1');
t2CondIdx = strcmp(behav.responseData_labels,'target type T2');
responseIdx = strcmp(behav.responseData_labels,'response');
correctIdx = strcmp(behav.responseData_labels,'correct');
rtIdx = strcmp(behav.responseData_labels,'RT');

cueCond = behav.responseData_all(:,cueCondIdx);
t1Cond = behav.responseData_all(:,t1CondIdx);
t2Cond = behav.responseData_all(:,t2CondIdx);
targetCond = [t1Cond t2Cond];
response = behav.responseData_all(:,responseIdx);
correct = behav.responseData_all(:,correctIdx);
rt = behav.responseData_all(:,rtIdx);

if any(strcmp(behav.responseData_labels, 'target pedestal T1'))
    pedestalOn = true;
    t1PedIdx = strcmp(behav.responseData_labels,'target pedestal T1');
    t2PedIdx = strcmp(behav.responseData_labels,'target pedestal T2');
    pedestal = behav.responseData_all(:,t1PedIdx | t2PedIdx);
else
    pedestalOn = false;
end

%% cue type (valid, invalid)
cueValidity = zeros(nTrials,1);
cueValidity(cueCond==2 | cueCond==5) = 1; % '1-1','2-2' valid
cueValidity(cueCond==3 | cueCond==4) = -1; % '1-2','2-1' invalid

%% target type (CCW, CW, absent)
responseTarget = zeros(nTrials,1);
responseTarget(cueCond==2 | cueCond==4) = 1; % '1-1','2-1'
responseTarget(cueCond==3 | cueCond==5) = 2; % '1-2','2-2'

targetType = nan(nTrials,1);
nontargetType = nan(nTrials,1);
targetPedestal = nan(nTrials,1);
nontargetPedestal = nan(nTrials,1);
for i = 1:nTrials
    if responseTarget(i)~=0
        targetType(i) = targetCond(i,responseTarget(i));
        nontargetType(i) = targetCond(i,3-responseTarget(i));
        if pedestalOn
            targetPedestal(i) = pedestal(i,responseTarget(i));
            nontargetPedestal(i) = pedestal(i,3-responseTarget(i));
        end
    end
end

% sanity check (compare to response and correct)
if pedestalOn
    correctResponse = (targetPedestal-1)*2 + targetType;
    swapResponse = (nontargetPedestal-1)*2 + nontargetType;
else
    correctResponse = targetType;
    correctResponse(targetType==0) = 3;
    swapResponse = nontargetType;
    swapResponse(nontargetType==0) = 3;
end

%% detection performance
targetPresent = nan(nTrials,1);
targetPresent(targetType==1 | targetType==2) = 1;
targetPresent(targetType==0) = 0;

% determine response options
responseOptions = unique(response(~isnan(response)));

presentResponse = nan(nTrials,1);
if numel(responseOptions)==4
    presentResponse(response==1 | response==2 | response==3 | response==4) = 1;
else
    presentResponse(response==1 | response==2) = 1;
    presentResponse(response==3) = 0;
end

detectHit = targetPresent==1 & presentResponse==1;
detectMiss = targetPresent==1 & presentResponse==0;
detectFA = targetPresent==0 & presentResponse==1;
detectCR = targetPresent==0 & presentResponse==0;

detectHMFC = double([detectHit detectMiss detectFA detectCR]);
detectHMFC(isnan(targetPresent),:) = NaN;
detectHMFC(targetPresent==0,1:2) = NaN; % target present trials only for H & M
detectHMFC(targetPresent==1,3:4) = NaN; % target absent trials only for FA & CR

%% discrimination performance
discrimCorrect = targetPresent==1 & correct==1;
discrimIncorrect = targetPresent==1 & presentResponse==1 & correct==-1;

discrimCI = double([discrimCorrect discrimIncorrect]);
discrimCI(isnan(targetPresent),:) = NaN;
discrimCI(targetPresent==0,:) = NaN;

%% overall accuracy
acc = behav.responseData_all(:,correctIdx);
acc(acc==-1) = 0;

%% exclude missed responses and wrong button presses
wWrongButton = behav.responseData_all(:,responseIdx)==0;
wNotMissed = behav.responseData_all(:,correctIdx)==1 | behav.responseData_all(:,correctIdx)==-1;
w = wWrongButton | ~wNotMissed;
detectHMFC(w,:) = NaN;
discrimCI(w,:) = NaN;
acc(w,:) = NaN;
rt(w,:) = NaN;

%% store
behav.cueValidity = cueValidity;
behav.responseTarget = responseTarget;
behav.targetType = targetType;
behav.nontargetType = nontargetType;
if pedestalOn
    behav.targetPedestal = targetPedestal;
    behav.nontargetPedestal = nontargetPedestal;
end
behav.correctResponse = correctResponse;
behav.swapResponse = swapResponse;
behav.response = response;
behav.targetPresent = targetPresent;
behav.presentResponse = presentResponse;
behav.detectHMFC = detectHMFC;
behav.discrimCI = discrimCI;
behav.acc = acc;
behav.rt = rt;

