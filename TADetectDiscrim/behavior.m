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

cueCond = behav.responseData_all(:,cueCondIdx);
t1Cond = behav.responseData_all(:,t1CondIdx);
t2Cond = behav.responseData_all(:,t2CondIdx);
targetCond = [t1Cond t2Cond];
response = behav.responseData_all(:,responseIdx);
correct = behav.responseData_all(:,correctIdx);

%% cue type (valid, invalid)
cueValidity = zeros(nTrials,1);
cueValidity(cueCond==2 | cueCond==5) = 1; % '1-1','2-2' valid
cueValidity(cueCond==3 | cueCond==4) = -1; % '1-2','2-1' invalid

%% target type (CCW, CW, absent)
responseTarget = zeros(nTrials,1);
responseTarget(cueCond==2 | cueCond==4) = 1; % '1-1','2-1'
responseTarget(cueCond==3 | cueCond==5) = 2; % '1-2','2-2'

targetType = nan(nTrials,1);
for i = 1:nTrials
    if responseTarget(i)~=0
        targetType(i) = targetCond(i,responseTarget(i));
    end
end

% sanity check (compare to response and correct)
correctResponse = targetType;
correctResponse(targetType==0) = 3;

%% detection performance
targetPresent = nan(nTrials,1);
targetPresent(targetType==1 | targetType==2) = 1;
targetPresent(targetType==0) = 0;

presentResponse = nan(nTrials,1);
presentResponse(response==1 | response==2) = 1;
presentResponse(response==3) = 0;

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

%% store
behav.cueValidity = cueValidity;
behav.responseTarget = responseTarget;
behav.targetType = targetType;
behav.targetPresent = targetPresent;
behav.presentResponse = presentResponse;
behav.detectHMFC = detectHMFC;
behav.discrimCI = discrimCI;
behav.acc = acc;

