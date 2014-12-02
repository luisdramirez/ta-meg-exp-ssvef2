function [Pc,responseData_all,responseData_labels] = TADetect_analysis
%% exp setup
trialCount = 41;  
respSecs = 1.4;
feedbackDur = 0.3;
refreshRate = 60;  %(frames)
blockLength = 300; %(frames)
respTime = 198;  % frames to respond period
keyCodes = [30 31];

%% combine responseData for all runs 
rootDir = '/Users/liusirui/Documents/MATLAB/MEG';
dataDir = [rootDir '/data/TADetect_data/TADetect_Rachel/'];
stimDir = [rootDir '/TA_EXP/TADetect/stimuli/TADetect_Rachel_Stim/'];
df = dir([dataDir,'*.mat']);
sf = dir([stimDir,'*.mat']);

if length(df) ~= length(sf)
    warning('inconsistent number of stim files and dataset')
end

responseData_all = [];

for n = 1:length(df); 
    dd = load([dataDir,df(n).name]);
    stim = load([stimDir,sf(n).name]);
    runNum = repmat(n,trialCount,1);
    [temp, responseData_labels] = sl_responseData(respTime,trialCount,...
        respSecs,feedbackDur,refreshRate, blockLength, keyCodes, dd.response, ...
        stim.order,n);
    responseData_all = [responseData_all;temp];
    
end

%% extract block order from responseData_all
run = responseData_all(:,1);
cueBlockOrder = responseData_all(:,4);
targetBlockOrder = responseData_all(:,6);
accuracy_all = responseData_all(:,9);

%% For each run calculate accuracy for: 
% postcue T1 valid & invalid; postcue T2 valid & invalid; overall accuracy

% blockNames = {'blank','fast-left'}; % fast-left
% attBlockNames = {'no-att','att-right'}; % att-right
% targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres','abs-abs'};
% cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'}; % 2-1 = cueT2,postcueT1

for n = 1:length(df)
    
    runIndx = run == n;
    
    % postCue = T1; cueBlockNames is valid '1-1' or invalid '2-1'

    valid_T1 = accuracy_all(runIndx & cueBlockOrder == 2);
    Pc.valid_T1(:,n) = sum (valid_T1 == 1) / numel(valid_T1);

    invalid_T1 = accuracy_all(runIndx & cueBlockOrder == 4);
    Pc.invalid_T1(:,n) = sum(invalid_T1 == 1) / numel(invalid_T1);

    % postCue = T2; cueBlockNames is valid '2-2 or invalid '1-2'

    valid_T2 = accuracy_all(runIndx & cueBlockOrder == 5);
    Pc.valid_T2(:,n) = sum(valid_T1 == 1) / numel(valid_T2);

    invalid_T2 = accuracy_all(runIndx & cueBlockOrder == 3);
    Pc.invalid_T2(:,n) = sum(invalid_T2 == 1) / numel(invalid_T2);

    % overall accuracy
    correctPerRun  = accuracy_all(runIndx) == 1;
    Pc.all(:,n) = sum(correctPerRun) / sum(responseData_all(runIndx, ...
        3) ~= 1);
end

%% calculate means and stes
% means
Pc.means = [mean(Pc.valid_T1),mean(Pc.invalid_T1),mean(Pc.valid_T2),...
    mean(Pc.invalid_T2)];

% stds
Pc.stds = [std(Pc.valid_T1),std(Pc.invalid_T1),std(Pc.valid_T2),...
    std(Pc.invalid_T2)];
% stes
Pc.stes = Pc.stds ./ sqrt (length(df));


%% plot
figure
hold on
% y = bar([Pc.means(1:2);Pc.means(3:4)],0.5);
y = errorbar([Pc.means(1:2);Pc.means(3:4)],[Pc.stes(1:2);Pc.stes(3:4)],'.');
ylim([0 1])
set(gca,'XTickLabel',{'','T1','','','','','T2',''});
ylabel('Accuracy')
legend(y,{'valid','invalid'});

end


