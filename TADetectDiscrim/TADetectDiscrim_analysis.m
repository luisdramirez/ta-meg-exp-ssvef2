function [Pc,responseData_all,responseData_labels] = TADetectDiscrim_analysis
%% exp setup
trialCount = 41;  
respSecs = 1.4;
feedbackDur = 0.3;
refreshRate = 60;  %(frames)
blockLength = 300; %(frames)
respTime = 198;  % frames to respond period
keyCodes = [30 31 32];

%% combine responseData for all runs 
rootDir = '/Users/liusirui/Documents/MATLAB/MEG/data/TADetectDiscrim/pilot';
dataDir = [rootDir '/rd/11-20/'];
stimDir = [rootDir '/rd/stimuli/11-20/'];
df = dir([dataDir,'*.mat']);
sf = dir([stimDir,'*.mat']);
%%
if length(df) ~= length(sf)
    warning('inconsistent number of stim files and dataset')
end

responseData_all = [];

for n = 1:length(df); 
    dd = load([dataDir,df(n).name]);
    stim = load([stimDir,sf(n).name]);
    runNum = repmat(n,trialCount,1);
    [temp, responseData_labels] = sl_responseDiscrimData(respTime,trialCount,...
        respSecs,feedbackDur,refreshRate, blockLength, keyCodes, dd.response, ...
        stim.order,n);
    responseData_all = [responseData_all;temp];
    
end

%% extract block order from responseData_all
run = responseData_all(:,1);
cueBlockOrder = responseData_all(:,4);
targetBlockOrder = responseData_all(:,6);
accuracy_all = responseData_all(:,9);
response_all = responseData_all(:,8);

%% For each run calculate detection and discrimination accuracy for: 
% postcue T1 valid & invalid; postcue T2 valid & invalid; overall accuracy

% blockNames = {'blank','fast-left'}; % fast-left
% attBlockNames = {'no-att','att-right'}; % att-right
% targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres','abs-abs'};
% cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'}; % 2-1 = cueT2,postcueT1

for n = 1:length(df)
    
    runIndx = run == n;
    
    % postCue = T1; cueBlockNames is valid '1-1' or invalid '2-1'

    valid_T1 = accuracy_all(runIndx & cueBlockOrder == 2);
    response_valid_T1 = response_all(runIndx & cueBlockOrder == 2);
    Pc.discrim_valid_T1(:,n) = sum (valid_T1 == 1 & ...
        response_valid_T1 ~= 3) / sum(response_valid_T1 ~= 3);
    Pc.detect_valid_T1(:,n) = sum (valid_T1 == 1) / numel(valid_T1);

    invalid_T1 = accuracy_all(runIndx & cueBlockOrder == 4);
    response_invalid_T1 = response_all(runIndx & cueBlockOrder == 4);
    Pc.discrim_invalid_T1(:,n) = sum (invalid_T1 == 1 & ...
        response_invalid_T1 ~= 3) / sum(response_invalid_T1 ~= 3);
    Pc.detect_invalid_T1(:,n) = sum(invalid_T1 == 1) / numel(invalid_T1);

    % postCue = T2; cueBlockNames is valid '2-2 or invalid '1-2'

    valid_T2 = accuracy_all(runIndx & cueBlockOrder == 5);
    response_valid_T2 = response_all(runIndx & cueBlockOrder == 5);
    Pc.discrim_valid_T2(:,n) = sum (valid_T2 == 1 & ...
        response_valid_T2 ~= 3) / sum(response_valid_T2 ~= 3);
    Pc.detect_valid_T2(:,n) = sum(valid_T2 == 1) / numel(valid_T2);

    invalid_T2 = accuracy_all(runIndx & cueBlockOrder == 3);
    response_invalid_T2 = response_all(runIndx & cueBlockOrder == 3);
    Pc.discrim_invalid_T2(:,n) = sum (invalid_T2 == 1 & ...
        response_invalid_T2 ~= 3) / sum(response_invalid_T2 ~= 3);
    Pc.detect_invalid_T2(:,n) = sum(invalid_T2 == 1) / numel(invalid_T2);

    % overall accuracy
    Detect_correctPerRun  = accuracy_all(runIndx) == 1;
    Discrim_correctPerRun = accuracy_all(runIndx) == 1 & ...
        response_all(runIndx) ~= 3;
    Pc.Detect_all(:,n) = sum(Detect_correctPerRun) /...
        sum(responseData_all(runIndx, 3) ~= 1);
    Pc.Discrim_all(:,n) = sum(Discrim_correctPerRun) / ...
        sum(responseData_all(runIndx,3) ~= 1 & response_all (runIndx)~= 3);
end

%% calculate means and stes
% means
Pc.discrim_means = [mean(Pc.discrim_valid_T1),mean(Pc.discrim_invalid_T1),mean(Pc.discrim_valid_T2),...
    mean(Pc.discrim_invalid_T2)];
Pc.detect_means = [mean(Pc.detect_valid_T1),mean(Pc.detect_invalid_T1),mean(Pc.detect_valid_T2),...
    mean(Pc.detect_invalid_T2)];

% stds
Pc.discrim_stds = [std(Pc.discrim_valid_T1),std(Pc.discrim_invalid_T1),std(Pc.discrim_valid_T2),...
    std(Pc.discrim_invalid_T2)];
Pc.detect_stds = [std(Pc.detect_valid_T1),std(Pc.detect_invalid_T1),std(Pc.detect_valid_T2),...
    std(Pc.detect_invalid_T2)];

% stes
Pc.discrim_stes = Pc.discrim_stds ./ sqrt (length(df));
Pc.detect_stes = Pc.detect_stds ./ sqrt (length(df));


%% plot
figure
hold on
subplot(2,1,1)
% y = bar([Pc.means(1:2);Pc.means(3:4)],0.5);
y = errorbar([Pc.discrim_means(1:2);Pc.discrim_means(3:4)],[Pc.discrim_stes(1:2);Pc.discrim_stes(3:4)],'.');
ylim([0 1])
set(gca,'XTickLabel',{'','T1','','','','','T2',''});
ylabel('Accuracy')
legend(y,{'valid','invalid'});
title('discrimination accuracy')

subplot(2,1,2)
y = errorbar([Pc.detect_means(1:2);Pc.detect_means(3:4)],[Pc.detect_stes(1:2);Pc.detect_stes(3:4)],'.');
ylim([0 1])
set(gca,'XTickLabel',{'','T1','','','','','T2',''});
ylabel('Accuracy')
legend(y,{'valid','invalid'});
title('detection accuracy')
end


