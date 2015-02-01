function [Pc,responseData_all,responseData_labels] = TADetectDiscrim_analysis(subject, runs, date)

% [Pc,responseData_all,responseData_labels] = TADetectDiscrim_analysis(subject, runs, [date])
% 
% subject is the subject ID, e.g. 'sl'
% runs is the run numbers to analyze, e.g. 1:10
% date (optional) is a string with the date of the experiment, e.g. '20141219'
%   if date is not given, will look for runs from all dates.

%% exp setup (now read from stim file, below)
% trialCount = 41;  
% respSecs = 1.4;
% feedbackDur = 0.3;
% refreshRate = 60;  %(frames)
% % blockLength = 300; %(frames) - 800ms SOA
% blockLength = (5-0.5)*60; %(frames) - other SOA (difference from 800)
% % respTime = 198;  % frames to respond period - 800ms SOA
% respTime = (5-1.7)*60; %frames to respond period - other SOA (difference from 800)
% keyCodes = [30 31 32];

%% add path
% addpath /Users/liusirui/Documents/MATLAB/MEG/ta-meg-exp/TADetectDiscrim

%% combine responseData for all runs 
% get the data from the server using pathToExpt
rootDir = pathToExpt;
% rootDir = pwd;
dataDir = sprintf('%s/data/%s', rootDir, subject);
stimDir = sprintf('%s/stimuli', rootDir);
% df = dir([dataDir,'*.mat']);

if ~exist('date','var')
    date = [];
end

for iRun = 1:numel(runs)
    run = runs(iRun);
    d = dir(sprintf('%s/%s*%d.mat', dataDir, date, run));
    if numel(d)~=1
        error('More or fewer than one matching data file: %s/*%d.mat', dataDir, run)
    else
        df(iRun) = d;
    end
end

%% load exp params from a sample run
% assumes all runs have identical params
expNum = regexp(df(1).name, '_taDetectDiscrim(\d*).mat','tokens');
stim = load(sprintf('%s/taDetectDiscrim%s.mat', stimDir, expNum{1}{1}));

trialCount = length(stim.p.blockOrder);  
respSecs = stim.p.respDur;
feedbackDur = stim.p.feedbackDur;
refreshRate = stim.p.refrate;  %(frames)
blockLength = stim.p.blockDur*refreshRate; %(frames)
respTime = blockLength-(respSecs+feedbackDur)*refreshRate; %frames to respond period
keyCodes = stim.p.keyCodes;

%% analyze each run
responseData_all = [];

for n = 1:length(df); 
    name = df(n).name;
    dd = load(sprintf('%s/%s',dataDir,name));
    expNum = regexp(name,'_taDetectDiscrim(\d*).mat','tokens');
    stim = load(sprintf('%s/taDetectDiscrim%s.mat', stimDir, expNum{1}{1}));
    [temp, responseData_labels] = sl_responseDiscrimData(respTime,trialCount,...
        respSecs,feedbackDur,refreshRate, blockLength, keyCodes, dd.response, ...
        stim.order,n);
    responseData_all = [responseData_all; temp];  
end

%% extract block order from responseData_all
run = responseData_all(:,1);
cueBlockOrder = responseData_all(:,4);
response_correct = responseData_all(:,11);
response_all = responseData_all(:,10);
targetTypeT1 = responseData_all(:,7);
targetTypeT2 = responseData_all(:,8);

% convert target type and response data for computing detection rate
DetectTargetType = responseData_all(:,7:8);
DetectTargetType (DetectTargetType == 2) = 1;
DetectResponse_all = response_all;
DetectResponse_all (DetectResponse_all  == 2) = 1;
DetectResponse_all (DetectResponse_all  == 3) = 0;


%% For each run calculate detection and discrimination accuracy for: 
% postcue T1 valid & invalid; postcue T2 valid & invalid; overall accuracy

% blockNames = {'blank','fast-left'}; % fast-left
% attBlockNames = {'no-att','att-right'}; % att-right
% targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres','abs-abs'};
% cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'}; % 2-1 = cueT2,postcueT1

for n = 1:length(df)
    
    runIndx = run == n;
    
    % postCue = T1; cueBlockNames is valid '1-1' or invalid '2-1'
    
    valid_T1_Indx = runIndx & cueBlockOrder == 2;   
    invalid_T1_Indx = runIndx & cueBlockOrder == 4;
    
     Pc.discrim_valid_T1(:,n) = sum (targetTypeT1(valid_T1_Indx) == response_all(valid_T1_Indx))/ ...
         sum(ismember(targetTypeT1(valid_T1_Indx),[1,2]) );
     Pc.detect_valid_T1(:,n) = sum (DetectTargetType(valid_T1_Indx,1) == DetectResponse_all(valid_T1_Indx) ) /...
         numel(DetectTargetType(valid_T1_Indx,1));
     Pc.overall_valid_T1(:,n) = sum ( response_correct(valid_T1_Indx) == 1 ) / ...
         numel( response_correct(valid_T1_Indx) );
     
     Pc.discrim_invalid_T1(:,n) = sum (targetTypeT1(invalid_T1_Indx) == response_all(invalid_T1_Indx))/ ...
         sum(ismember(targetTypeT1(invalid_T1_Indx),[1,2]) );
     Pc.detect_invalid_T1(:,n) = sum (DetectTargetType(invalid_T1_Indx,1) == DetectResponse_all(invalid_T1_Indx) ) /...
         numel(DetectTargetType(invalid_T1_Indx,1));
     Pc.overall_invalid_T1(:,n) = sum ( response_correct(invalid_T1_Indx) == 1) / ...
         numel( response_correct(invalid_T1_Indx) );

    % postCue = T2; cueBlockNames is valid '2-2 or invalid '1-2'
     valid_T2_Indx = runIndx & cueBlockOrder == 5;
     invalid_T2_Indx = runIndx & cueBlockOrder == 3;

     Pc.discrim_valid_T2(:,n) = sum (targetTypeT2(valid_T2_Indx) == response_all(valid_T2_Indx))/ ...
         sum(ismember(targetTypeT2(valid_T2_Indx),[1,2]) );
     Pc.detect_valid_T2(:,n) = sum (DetectTargetType(valid_T2_Indx,2) == DetectResponse_all(valid_T2_Indx) ) /...
         numel(DetectTargetType(valid_T2_Indx,2));
     Pc.overall_valid_T2(:,n) = sum( response_correct(valid_T2_Indx) == 1 ) /...
         numel( response_correct(valid_T2_Indx) );

     Pc.discrim_invalid_T2(:,n) = sum( targetTypeT2(invalid_T2_Indx) == response_all(invalid_T2_Indx) )/ ...
         sum( ismember(targetTypeT2(invalid_T2_Indx),[1,2]) );
     Pc.detect_invalid_T2(:,n) = sum (DetectTargetType(invalid_T2_Indx,2) == DetectResponse_all(invalid_T2_Indx) ) /...
         numel(DetectTargetType(invalid_T2_Indx,2));
     Pc.overall_invalid_T2(:,n) = sum( response_correct(invalid_T2_Indx) == 1) /...
         numel( response_correct(invalid_T2_Indx) );

     
     
    %% overall accuracy (detect and discrim seperate)
%     cueT1 = runIndx & (valid_T1_Indx | invalid_T1_Indx);
%     cueT2 = runIndx & (valid_T2_Indx | invalid_T2_Indx);
%     
%     pc.DetectT1(:,n) = sum( DetectTargetType (cueT1,1) == ...
%         DetectResponse_all (cueT1) ) / numel(DetectTargetType (cueT1,1));
%     
%     pc.DetectT2(:,n) = sum( DetectTargetType (cueT2,2) == ...
%         DetectResponse_all (cueT2) ) / numel (DetectTargetType (cueT2,2));
% 
%     pc.DiscrimT1(:,n) = sum( targetTypeT1 (cueT1) == ...
%         response_all(cueT1) ) / sum ( ismember(targetTypeT1(cueT1),[1,2]) );
%        
%     pc.DiscrimT2(:,n) = sum( targetTypeT2 (cueT2) == ...
%         response_all(cueT2) ) / sum( ismember(targetTypeT2(cueT2),[1,2]) );
  
end



%% calculate means and stes
% means
Pc.discrim_means = [mean(Pc.discrim_valid_T1),mean(Pc.discrim_invalid_T1),mean(Pc.discrim_valid_T2),...
    mean(Pc.discrim_invalid_T2)];
Pc.detect_means = [mean(Pc.detect_valid_T1),mean(Pc.detect_invalid_T1),mean(Pc.detect_valid_T2),...
    mean(Pc.detect_invalid_T2)];
Pc.overall_mean = [mean(Pc.overall_valid_T1), mean(Pc.overall_invalid_T1), mean(Pc.overall_valid_T2),...
    mean(Pc.overall_invalid_T2) ];

% Pc.Discrim_all_mean = [mean(pc.DiscrimT1),mean(pc.DiscrimT2)];
% Pc.Detect_all_mean = [mean(pc.DetectT1),mean(pc.DetectT2)];

% stds
Pc.discrim_stds = [std(Pc.discrim_valid_T1),std(Pc.discrim_invalid_T1),std(Pc.discrim_valid_T2),...
    std(Pc.discrim_invalid_T2)];
Pc.detect_stds = [std(Pc.detect_valid_T1),std(Pc.detect_invalid_T1),std(Pc.detect_valid_T2),...
    std(Pc.detect_invalid_T2)];
Pc.overall_stds = [std(Pc.overall_valid_T1), std(Pc.overall_invalid_T1), std(Pc.overall_valid_T2),...
    std(Pc.overall_invalid_T2) ];

% Pc.Discrim_all_std = [std(pc.DiscrimT1),std(pc.DiscrimT2)];
% Pc.Detect_all_std = [std(pc.DetectT1),std(pc.DetectT2)];

% stes
Pc.discrim_stes = Pc.discrim_stds ./ sqrt (length(df));
Pc.detect_stes = Pc.detect_stds ./ sqrt (length(df));
Pc.overall_stes = Pc.overall_stds ./ sqrt( length(df) );

% Pc.Discrim_all_ste = Pc.Discrim_all_std ./ sqrt (length(df));
% Pc.Detect_all_ste = Pc.Detect_all_std ./ sqrt(length(df));

%% plot
scrsz=get(0,'ScreenSize');                                 
% figure('Position', [1 scrsz(4)*1/3 scrsz(3)*1/3 scrsz(4)])
figure

hold on
subplot(1,3,1)
% y = bar([Pc.means(1:2);Pc.means(3:4)],0.5);
y = errorbar([Pc.discrim_means(1:2);Pc.discrim_means(3:4)],[Pc.discrim_stes(1:2);Pc.discrim_stes(3:4)],'.');
ylim([0 1])
% set(y, 'MarkerSize', 20)
set(y(2),'Color','r')
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'T1','T2'});
% set(gca,'XTickLabel',{'','T1','','','','','T2',''});
ylabel('accuracy')
% legend(y,{'valid','invalid'});
legend('valid','invalid','Location','SouthEast');
title('discrimination')

subplot(1,3,2)
y = errorbar([Pc.detect_means(1:2);Pc.detect_means(3:4)],[Pc.detect_stes(1:2);Pc.detect_stes(3:4)],'.');
ylim([0 1])
% set(y, 'MarkerSize', 20)
set(y(2),'Color','r')
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'T1','T2'});
% set(gca,'XTickLabel',{'','T1','','','','','T2',''});
ylabel('accuracy')
title('detection')

subplot(1,3,3)
y = errorbar([Pc.overall_mean(1:2);Pc.overall_mean(3:4)],[Pc.overall_stes(1:2);Pc.overall_stes(3:4)],'.');
ylim([0 1])
% set(y, 'MarkerSize', 20)
set(y(2),'Color','r')
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'T1','T2'});
% set(gca,'XTickLabel',{'','T1','','','','','T2',''});
ylabel('accuracy')
title('overall')

% turnwhite

%%
% figure
% y = errorbar([Pc.Discrim_all_mean;Pc.Detect_all_mean],[Pc.Discrim_all_ste;Pc.Detect_all_ste],'.');
% ylim([0 1])
% set(gca,'XTickLabel',{'','discrimination','','','','','detection',''});
% legend(y,{'T1','T2'},'Location','NorthWest')
% ylabel('Accuracy')
% title('overall accuracy')

end


