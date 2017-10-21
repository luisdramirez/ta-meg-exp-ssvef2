function [accuracy,accuracy2,responseData_all,responseData_labels, io] = TADetectDiscrim_analysis(subject, runs, dates, plotLevel,saveFile,saveFigs)

% function [accuracy,accuracy2,responseData_all,responseData_labels] = TADetectDiscrim_analysis(subject, runs, dates, plotLevel,saveFile,saveFigs)
%
% subject is the subject ID, e.g. 'sl'
% runs is the run numbers to analyze, e.g. 1:10
% date (optional) is a string with the date of the experiment, e.g. '20141219'
%   if date is not given, will look for runs from all dates.

if nargin < 3 || isempty(dates)
    dates = {[]};
end

if nargin < 4 || isempty(plotLevel)
    plotLevel = 1; % 1 plot minimal --> 3 plot everything
end

if nargin < 5 || isempty(saveFile)
    saveFile = 0;
end

if nargin < 6 || isempty(saveFigs)
    saveFigs = 0;
end

if ischar(dates)
    dateStr = dates;
    clear dates
    dates{1} = dateStr;
end

subjectStr = sprintf('%s_taDetectDiscrim', subject);
% subjectStr = 'taDetectDiscrim';

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
rootDir = pathToMEGExpt;
% rootDir = pathToExpt;
% rootDir = '~/Desktop/Luis/ta-meg-exp-ssvef2/TADetectDiscrim';
% rootDir = '/Volumes/purplab/EXPERIMENTS/1_Current Experiments/Luis/ta-meg-exp-ssvef2/TADetectDiscrim';
%rootDir = '/Users/luisramirez/Documents/CarrascoLabMEG/ta-meg-exp-ssvef2/TADetectDiscrim';
% rootDir = pwd;
dataDir = sprintf('%s/data/%s', rootDir, subject);
stimDir = sprintf('%s/stimuli', rootDir);
analysisDir = sprintf('%s/analysis/%s',rootDir,subject);
figDir = sprintf('%s/figures/%s',rootDir,subject);

if ~exist(analysisDir,'dir')
    mkdir(analysisDir);
end
if ~exist(figDir,'dir')
    mkdir(figDir)
end

if numel(dates)==1
    dateStr = dates{1};
else
    dateStr = sprintf('%s-%s', dates{1}, dates{end});
end
% df = dir([dataDir,'*.mat']);
if numel(runs)>1
    analysisFile = sprintf('%s_%s_taDetectDiscrim_%s%s%s',subject,dateStr,num2str(runs(1)),'-',num2str(runs(end)));
else
    analysisFile = sprintf('%s_%s_taDetectDiscrim_%s',subject,dateStr,num2str(runs));
end

counter = 1;
for iDate = 1:numel(dates)
    date = dates{iDate};
    for iRun = 1:numel(runs)
        run = runs(iRun);
        d = dir(sprintf('%s/%s*taDetectDiscrim%d.mat', dataDir, date, run));
        if numel(d)~=1
            error('More or fewer than one matching data file: %s/*%d.mat', dataDir, run)
        else
            %             df(iRun) = d;
            df(counter) = d;
            counter = counter + 1;
        end
    end
end

%% load exp params from a sample run
% assumes all runs have identical params
expNum = regexp(df(1).name, '_taDetectDiscrim(\d*).mat','tokens');
stim = load(sprintf('%s/%s%s.mat', stimDir, subjectStr, expNum{1}{1}));

trialCount = length(stim.p.blockOrder);
respSecs = stim.p.respDur;
feedbackDur = stim.p.feedbackDur;
refreshRate = stim.p.refrate;  %(frames)
blockLength = stim.p.blockDur*refreshRate; %(frames)
% respTime = blockLength-(respSecs+feedbackDur)*refreshRate; %frames to respond period
respTime = (stim.p.targetLeadTime + stim.p.targetSOA + stim.p.cueTargetSOA)*refreshRate;
keyCodes = stim.p.keyCodes;
responseOption = stim.p.responseOption;

%% analyze each run
responseData_all = [];
targetPosType_all = [];

for n = 1:length(df)
    name = df(n).name;
    dd = load(sprintf('%s/%s',dataDir,name));
    expNum = regexp(name,'_taDetectDiscrim(\d*).mat','tokens');
    stim = load(sprintf('%s/%s%s.mat', stimDir, subjectStr, expNum{1}{1}));
    if isfield(stim.stimulus,'itiSeq')
        itiSeq = stim.stimulus.itiSeq;
    else
        itiSeq = [];
    end
    if isfield(stim.stimulus,'jitSeq')
        jitSeq = stim.stimulus.jitSeq;
    else
        jitSeq = [];
    end
    
    [temp, responseData_labels] = sl_responseDiscrimData(respTime,trialCount,...
        respSecs,refreshRate, blockLength, keyCodes, dd.response, ...
        stim.order, n, itiSeq, jitSeq, stim.p.jitter);
    responseData_all = [responseData_all; temp];
    
    % target pos type
    targetPos = stim.order.posBlockOrder';
    targetPosType = targetPos;
    targetPosType(targetPos>=1 & targetPos<=4) = 1;
    targetPosType(targetPos>=5 & targetPos<=8) = 2;
    targetPosType_all = [targetPosType_all; targetPosType];
    
    % target contrast
    if isfield(dd.response.target,'contrast')
        targetContrast(n,:) = dd.response.target.contrast;
    end
end

%% extract block order from responseData_all
run = responseData_all(:,1);
cueBlockOrder = responseData_all(:,4);
targetCondition = responseData_all(:,6);
targetTypeT1 = responseData_all(:,7);
targetTypeT2 = responseData_all(:,8);
response_all = responseData_all(:,10);
response_correct = responseData_all(:,11);
response_rt = responseData_all(:,12);
targetAxis = responseData_all(:,13:14);
targetPedestal = responseData_all(:,15:16);

% convert target type and response data for computing detection rate
DetectTargetType = responseData_all(:,7:8);
DetectTargetType (DetectTargetType == 2) = 1;
DetectResponse_all = response_all;
DiscrimResponse_all = response_all;
OverallResponse_all = response_all;
        
switch stim.p.responseOption
    case 'targetContrast4Levels'
        DetectResponse_all (DetectResponse_all >= 1 & DetectResponse_all <= 4) = 1;
        DiscrimResponse_all (DiscrimResponse_all >= 1 & DiscrimResponse_all <= 4) = 1;
    otherwise
        DetectResponse_all (DetectResponse_all  == 2) = 1;
        DetectResponse_all (DetectResponse_all  == 3) = 0;
        DiscrimResponse_all (DiscrimResponse_all  == 2) = 1;
        OverallResponse_all ( OverallResponse_all == 3) = 0;
end

%% For each run calculate detection and discrimination accuracy for:
% postcue T1 valid & invalid; postcue T2 valid & invalid; overall accuracy

% blockNames = {'blank','fast-left'}; % fast-left
% attBlockNames = {'no-att','att-right'}; % att-right
% targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres','abs-abs'};
% cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'}; % 2-1 = cueT2,postcueT1

cueBlockOrder_Indx = [2 4 5 3];
accuracy =[]; % rows: '1-1','2-1','2-2','1-2'; columns: runs 1-10
accuracy.missedTrials = zeros(1,length(df)); % columns: runs 1-10

for n = 1:length(df)
    runIndx = run == n;
    for k = 1:length(cueBlockOrder_Indx) % cueBlockOrder_Indx = [2 4 5 3]; cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'};
        condition_Indx = runIndx & cueBlockOrder == cueBlockOrder_Indx(k);
        if cueBlockOrder_Indx(k) == 2 || cueBlockOrder_Indx(k) == 4
            interval = 1;
            switch responseOption
                case 'targetType'
                    targetType = targetTypeT1;
                case 'targetContrast4Levels'
                    if numel(stim.p.keyCodes==4)
                        targetType = 2*(targetPedestal(:,1)-1) + targetTypeT1;
                    else
                        targetType = targetTypeT1;
                    end
                case 'targetPos'
                    targetType = targetPosType_all(:,1);
                otherwise
                    error('responseOption not recognized')
            end
        else
            interval = 2;
            switch responseOption
                case 'targetType'
                    targetType = targetTypeT2;
                case 'targetContrast4Levels'
                    if numel(stim.p.keyCodes==4)
                        targetType = 2*(targetPedestal(:,2)-1) + targetTypeT2;
                    else
                        targetType = targetTypeT2;
                    end
                case 'targetPos'
                    targetType = targetPosType_all(:,2);
                otherwise
                    error('responseOption not recognized')
            end
        end
        % total correct / total present
        accuracy.Discrim_all(k,n) = sum (targetType(condition_Indx) == response_all(condition_Indx))/ ...
            sum(ismember(targetType(condition_Indx),[1,2]) );
        % total correct / total correctly detected
        accuracy.Discrim1_all(k,n) = sum (targetType(condition_Indx) == response_all(condition_Indx))/ ...
            sum (DetectTargetType(condition_Indx,interval) == DiscrimResponse_all(condition_Indx) );
        accuracy.Detect_all(k,n) = sum (DetectTargetType(condition_Indx,interval) == DetectResponse_all(condition_Indx) ) /...
            numel(DetectTargetType(condition_Indx,interval));
        accuracy.Overall_all(k,n) = sum ( response_correct(condition_Indx) == 1 ) / ...
            numel( response_correct(condition_Indx) );
        accuracy.Hit_all(k,n) = sum (DetectTargetType(condition_Indx,interval) == DiscrimResponse_all(condition_Indx) ) / ...
            sum(ismember(targetType(condition_Indx),[1,2]) );
        accuracy.FA_all(k,n) = sum (DetectTargetType(condition_Indx,interval)== 0 & DetectResponse_all(condition_Indx) == 1) / ...
            sum(ismember(targetType(condition_Indx),[1,2]) );
        accuracy.Miss_all(k,n) = 1 - accuracy.Hit_all(k,n);
        accuracy.CR_all(k,n) = 1 -  accuracy.FA_all(k,n);
        
        accuracy.Hit_alltemp(k,n) = accuracy.Hit_all(k,n);
        accuracy.FA_alltemp(k,n) = accuracy.FA_all(k,n);
        if accuracy.Hit_alltemp(k,n) == 1
            accuracy.Hit_alltemp(k,n) = .99; % avoid Inf for dprime
        elseif accuracy.Hit_alltemp(k,n) == 0
            accuracy.Hit_alltemp(k,n) = .01;
        end
        if accuracy.FA_alltemp(k,n) == 1
            accuracy.FA_alltemp(k,n) = .99;
        elseif accuracy.FA_alltemp(k,n) == 0
            accuracy.FA_alltemp(k,n) = .01;
        end
        accuracy.dprime(k,n) = norminv(accuracy.Hit_alltemp(k,n)) - norminv(accuracy.FA_alltemp(k,n));
        
        accuracy.RT(k,n) = nanmean(response_rt(condition_Indx));
        
        accuracy.missedTrials(n) = accuracy.missedTrials(n) + sum(isnan(response_all(condition_Indx)) );
    end
end


%% calculate means and stes
accuracy.Discrim_means = nanmean(accuracy.Discrim_all,2); % rows: '1-1','2-1','2-2','1-2'
accuracy.Discrim1_means = nanmean(accuracy.Discrim1_all,2);
accuracy.Detect_means = nanmean(accuracy.Detect_all,2);
accuracy.Overall_means = nanmean(accuracy.Overall_all,2);
accuracy.Hit_means = nanmean(accuracy.Hit_all,2);
accuracy.FA_means = nanmean(accuracy.FA_all,2);
accuracy.Miss_means = nanmean(accuracy.Miss_all,2);
accuracy.CR_means = nanmean(accuracy.CR_all,2);
accuracy.dprime_means = nanmean(accuracy.dprime,2);
accuracy.RT_means = nanmean(accuracy.RT,2);
all_means = [accuracy.Hit_means,accuracy.FA_means,accuracy.Miss_means,accuracy.CR_means ]; % rows: '1-1','2-1','2-2','1-2'
% columns: Hit,FA, Miss, CR
accuracy.Discrim_stes = nanstd(accuracy.Discrim_all,0,2)./sqrt(length(df));
accuracy.Discrim1_stes = nanstd(accuracy.Discrim1_all,0,2)./sqrt(length(df));
accuracy.Detect_stes = nanstd(accuracy.Detect_all,0,2)./sqrt(length(df));
accuracy.Overall_stes = nanstd(accuracy.Overall_all,0,2) ./ sqrt(length(df));
accuracy.Hit_stes = nanstd(accuracy.Hit_all,0,2) ./ sqrt(length(df));
accuracy.FA_stes = nanstd(accuracy.FA_all,0,2)./ sqrt(length(df));
accuracy.Miss_stes = nanstd(accuracy.Miss_all,0,2)./ sqrt(length(df));
accuracy.CR_stes = nanstd(accuracy.CR_all,0,2)./ sqrt(length(df));
accuracy.dprime_stes = nanstd(accuracy.dprime,0,2) ./ sqrt(length(df));
accuracy.RT_stes = nanstd(accuracy.RT,0,2) ./ sqrt(length(df));
all_stes = [accuracy.Hit_stes,accuracy.FA_stes,accuracy.Miss_stes,accuracy.CR_stes];

if exist('targetContrast','var')
    accuracy.targetContrast = targetContrast;
end

%% plot
switch plotLevel
    case {1,2}
        scrsz=get(0,'ScreenSize');
        f(1) = figure('Position', [1 scrsz(4) scrsz(3)*(3/5) scrsz(4)/2]);
        
        subplot(1,4,1)
        y = errorbar([ (accuracy.Detect_means(1:2)');(accuracy.Detect_means(3:4)')],[(accuracy.Detect_stes(1:2)');(accuracy.Detect_stes(3:4)')],'.');
        xlim([0.5 2.5])
        ylim([0 1])
        set(y,'MarkerSize',25)
        set(y(2),'Color','r')
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'});
        % set(gca,'XTickLabel',{'','T1','','','','','T2',''});
        ylabel('accuracy')
        title('detection')
        
        subplot(1,4,2)
        y = errorbar([ (accuracy.Discrim1_means(1:2)'); (accuracy.Discrim1_means(3:4)')],[(accuracy.Discrim1_stes(1:2)');(accuracy.Discrim1_stes(3:4)')],'.');
        xlim([0.5 2.5])
        ylim([0 1])
        set(y,'MarkerSize',25)
        set(y(2),'Color','r')
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'});
        ylabel('accuracy')
        title('discrimination(total correct/total detected)')
        
        subplot(1,4,3)
        y = errorbar([ (accuracy.Overall_means(1:2)');(accuracy.Overall_means(3:4)')],[(accuracy.Overall_stes(1:2)');(accuracy.Overall_stes(3:4)')],'.');
        xlim([0.5 2.5])
        ylim([0 1])
        set(y,'MarkerSize',25)
        set(y(2),'Color','r')
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'});
        ylabel('accuracy')
        title('overall')
        
        subplot(1,4,4)
        y = errorbar([ (accuracy.RT_means(1:2)');(accuracy.RT_means(3:4)')],[(accuracy.RT_stes(1:2)');(accuracy.RT_stes(3:4)')],'.');
        xlim([0.5 2.5])
        ylim([0 1.6])
        set(y,'MarkerSize',25)
        set(y(2),'Color','r')
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'});
        ylabel('RT')
        title('overall')
        legend('valid','invalid')
        
    case 3
        scrsz=get(0,'ScreenSize');
        f(1) = figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)/2]);
        
        hold on
        subplot(1,5,1)
        % y = bar([Pc.means(1:2);Pc.means(3:4)],0.5);
        y = errorbar([ (accuracy.Discrim_means (1:2)') ; (accuracy.Discrim_means (3:4)')],[(accuracy.Discrim_stes(1:2)');(accuracy.Discrim_stes(3:4)')],'.');
        ylim([0 1])
        % set(y, 'MarkerSize', 20)
        set(y(2),'Color','r')
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'});
        % set(gca,'XTickLabel',{'','T1','','','','','T2',''});
        ylabel('accuracy')
        % legend(y,{'valid','invalid'});
        legend('valid','invalid','Location','NorthEast');
        title('discrimination(total correct/total present)')
        
        subplot(1,5,2)
        % y = bar([Pc.means(1:2);Pc.means(3:4)],0.5);
        y = errorbar([ (accuracy.Discrim1_means(1:2)'); (accuracy.Discrim1_means(3:4)')],[(accuracy.Discrim1_stes(1:2)');(accuracy.Discrim1_stes(3:4)')],'.');
        ylim([0 1])
        % set(y, 'MarkerSize', 20)
        set(y(2),'Color','r')
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'});
        % set(gca,'XTickLabel',{'','T1','','','','','T2',''});
        ylabel('accuracy')
        % legend(y,{'valid','invalid'});
        title('discrimination(total correct/total detected)')
        
        subplot(1,5,3)
        y = errorbar([ (accuracy.Detect_means(1:2)');(accuracy.Detect_means(3:4)')],[(accuracy.Detect_stes(1:2)');(accuracy.Detect_stes(3:4)')],'.');
        ylim([0 1])
        % set(y, 'MarkerSize', 20)
        set(y(2),'Color','r')
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'});
        % set(gca,'XTickLabel',{'','T1','','','','','T2',''});
        ylabel('accuracy')
        title('detection')
        
        subplot(1,5,4)
        y = errorbar([ (accuracy.dprime_means(1:2)');(accuracy.dprime_means(3:4)') ],[(accuracy.dprime_stes(1:2)'); (accuracy.dprime_stes(3:4)')],'.');
        ylim([0 5])
        set(y(2),'Color','r')
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'});
        ylabel('d''')
        title('detection')
        
        subplot(1,5,5)
        y = errorbar([ (accuracy.Overall_means(1:2)');(accuracy.Overall_means(3:4)')],[(accuracy.Overall_stes(1:2)');(accuracy.Overall_stes(3:4)')],'.');
        ylim([0 1])
        % set(y, 'MarkerSize', 20)
        set(y(2),'Color','r')
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'});
        % set(gca,'XTickLabel',{'','T1','','','','','T2',''});
        ylabel('accuracy')
        title('overall')
    otherwise
        %         error('plotLevel not recognized')
end


%% plot (hit, fa, miss, cr)
if plotLevel>2
    f(2) = figure('Position', [1 1 scrsz(3)*3/4 scrsz(4)/2]);
    subplot(1,2,1)
    barwitherr ([all_stes(1,:)' all_stes(2,:)'],[1 2 3 4],[all_means(1,:)' all_means(2,:)'])
    ylim([0 1])
    legend('valid','invalid','Location','NorthEast')
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1]
    colormap(barmap);
    title('T1 detection')
    set(gca, 'XTick',[1 2 3 4],'XTickLabel',{'Hit','FA','Miss','CR' });
    
    subplot(1,2,2)
    barwitherr ([all_stes(3,:)' all_stes(4,:)'],[1 2 3 4],[all_means(3,:)' all_means(4,:)'])
    ylim([0 1])
    legend('valid','invalid','Location','NorthEast')
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1];
    colormap(barmap);
    title('T2 detection')
    set(gca, 'XTick',[1 2 3 4],'XTickLabel',{'Hit','FA','Miss','CR' });
end

%% pp pa ap aa
accuracy2 = []; % rows: pp,pa,ap,aa; columns: '1-1','2-1','2-2','1-2'; pages: runs 1-10

for n = 1:length(df)
    runIndx = run == n;
    for k = 1:length(cueBlockOrder_Indx) % cueBlockOrder_Indx = [2 4 5 3]; cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'};
        dd = runIndx & cueBlockOrder == cueBlockOrder_Indx(k);
        for i = 2:5 % targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres','abs-abs'};
            condition_Indx =  dd & targetCondition == i ;
            if cueBlockOrder_Indx(k) == 2 || cueBlockOrder_Indx(k) == 4
                interval = 1;
                targetType = targetTypeT1;
            else
                interval = 2;
                targetType = targetTypeT2;
            end
            % total correct / total present
            accuracy2.Discrim_all(i-1,k,n) = sum (targetType(condition_Indx) == response_all(condition_Indx))/ ...
                sum(ismember(targetType(condition_Indx),[1,2]) );
            % total correct / total correctly detected
            accuracy2.Discrim1_all(i-1,k,n) = sum (targetType(condition_Indx) == response_all(condition_Indx))/ ...
                sum (DetectTargetType(condition_Indx,interval) == DiscrimResponse_all(condition_Indx) );
            accuracy2.Detect_all(i-1,k,n) = sum (DetectTargetType(condition_Indx,interval) == DetectResponse_all(condition_Indx) ) /...
                numel(DetectTargetType(condition_Indx,interval));
            accuracy2.Overall_all(i-1,k,n) = sum ( response_correct(condition_Indx) == 1 ) / ...
                numel( response_correct(condition_Indx) );
            
        end
    end
end



%% calculate means and stes (pp pa ap aa)
accuracy2.Discrim_means = nanmean(accuracy2.Discrim_all,3); % rows:  pp,pa,ap,aa; columns: target type
accuracy2.Discrim1_means = nanmean(accuracy2.Discrim1_all,3);
accuracy2.Detect_means = nanmean(accuracy2.Detect_all,3);
accuracy2.Overall_means = nanmean(accuracy2.Overall_all,3);

accuracy2.Discrim_stes = nanstd(accuracy2.Discrim_all,0,3)./sqrt(length(df));
accuracy2.Discrim1_stes = nanstd(accuracy2.Discrim1_all,0,3)./sqrt(length(df));
accuracy2.Detect_stes = nanstd(accuracy2.Detect_all,0,3)./sqrt(length(df));
accuracy2.Overall_stes = nanstd(accuracy2.Overall_all,0,3) ./ sqrt(length(df));

accuracy2.T1_Discrim_means = accuracy2.Discrim_means(1:2,1:2); % rows: pp,pa; columns: T1 valid, invalid;
accuracy2.T1_Discrim_stes = accuracy2.Discrim_stes(1:2,1:2);
accuracy2.T2_Discrim_means = [accuracy2.Discrim_means(1,3:4) ; accuracy2.Discrim_means(3,3:4)]; % rows: pp,ap columns: T2 valid, invalid;
accuracy2.T2_Discrim_stes = [accuracy2.Discrim_stes(1,3:4) ; accuracy2.Discrim_stes(3,3:4)];

accuracy2.T1_Discrim1_means = accuracy2.Discrim1_means(1:2,1:2); % rows: pp,pa; columns:  T1 valid, invalid;
accuracy2.T1_Discrim1_stes = accuracy2.Discrim1_stes(1:2,1:2);
accuracy2.T2_Discrim1_means = [accuracy2.Discrim1_means(1,3:4) ; accuracy2.Discrim1_means(3,3:4)];% rows: pp,ap columns: T2 valid, invalid;
accuracy2.T2_Discrim1_stes = [accuracy2.Discrim1_stes(1,3:4) ; accuracy2.Discrim1_stes(3,3:4)];

accuracy2.T1_Detect_means = accuracy2.Detect_means(:,1:2);
accuracy2.T1_Detect_stes = accuracy2.Detect_stes(:,1:2);
accuracy2.T2_Detect_means = accuracy2.Detect_means(:,3:4);
accuracy2.T2_Detect_stes = accuracy2.Detect_stes(:,3:4);

accuracy2.T1_Overall_means = accuracy2.Overall_means(:,1:2);
accuracy2.T1_Overall_stes = accuracy2.Overall_stes(:,1:2);
accuracy2.T2_Overall_means = accuracy2.Overall_means(:,3:4);
accuracy2.T2_Overall_stes = accuracy2.Overall_stes(:,3:4);


%% plot discrimination pp pa/pp ap for T1 and T2
if plotLevel > 1
    f(4) = figure('Position', [1 scrsz(1) scrsz(3)*3/4 scrsz(4)/2]);
    subplot(1,2,1)
    
    barwitherr ([accuracy2.T1_Discrim1_stes(:,1) accuracy2.T1_Discrim1_stes(:,2)],[1 2],[accuracy2.T1_Discrim1_means(:,1) accuracy2.T1_Discrim1_means(:,2)])
    legend('valid','invalid','Location','NorthEast')
    ylim([0 1])
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1]
    colormap(barmap);
    title('T1 Discrim (total correct/total detected)')
    set(gca, 'XTick',[1 2],'XTickLabel',{'pp','pa' });
    
    subplot(1,2,2)
    barwitherr ([accuracy2.T2_Discrim1_stes(:,1) accuracy2.T2_Discrim1_stes(:,2)],[1 2],[accuracy2.T2_Discrim1_means(:,1) accuracy2.T2_Discrim1_means(:,2)])
    legend('valid','invalid','Location','NorthEast')
    ylim([0 1])
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1];
    colormap(barmap);
    title('T2 Discrim (total correct/total detected)')
    set(gca, 'XTick',[1 2],'XTickLabel',{'pp','ap' });
end

if plotLevel > 2
    f(3) = figure('Position', [1 scrsz(4) scrsz(3)*3/4 scrsz(4)/2]);
    subplot(1,2,1)
    ylim([0 1])
    barwitherr ([accuracy2.T1_Discrim_stes(:,1) accuracy2.T1_Discrim_stes(:,2)],[1 2],[accuracy2.T1_Discrim_means(:,1) accuracy2.T1_Discrim_means(:,2)])
    legend('valid','invalid','Location','best')
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1]
    colormap(barmap);
    title('T1 Discrim (total correct/total present)')
    set(gca, 'XTick',[1 2],'XTickLabel',{'pp','pa' });
    
    subplot(1,2,2)
    barwitherr ([accuracy2.T2_Discrim_stes(:,1) accuracy2.T2_Discrim_stes(:,2)],[1 2],[accuracy2.T2_Discrim_means(:,1) accuracy2.T2_Discrim_means(:,2)])
    legend('valid','invalid','Location','best')
    ylim([0 1])
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1];
    colormap(barmap);
    title('T2 Discrim (total correct/total present)')
    set(gca, 'XTick',[1 2],'XTickLabel',{'pp','ap' });
end

%% plot detection for pp ap ap aa
if plotLevel > 1
    f(5) = figure('Position', [1 1 scrsz(3)*3/4 scrsz(4)/2]);
    subplot(1,2,1)
    barwitherr ([accuracy2.T1_Detect_stes(:,1) accuracy2.T1_Detect_stes(:,2)],[1 2 3 4],[accuracy2.T1_Detect_means(:,1) accuracy2.T1_Detect_means(:,2)])
    legend('valid','invalid')
    ylim([0 1])
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1]
    colormap(barmap);
    title('T1 detection')
    set(gca, 'XTick',[1 2 3 4],'XTickLabel',{'pp','pa','ap','aa' });
    
    subplot(1,2,2)
    barwitherr ([accuracy2.T2_Detect_stes(:,1) accuracy2.T2_Detect_stes(:,2)],[1 2 3 4],[accuracy2.T2_Detect_means(:,1) accuracy2.T2_Detect_means(:,2)])
    legend('valid','invalid')
    ylim([0 1])
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1];
    colormap(barmap);
    title('T2 detection')
    set(gca, 'XTick',[1 2 3 4],'XTickLabel',{'pp','pa','ap','aa' });
end

%% plot overall accuracy for pp pa ap aa
if plotLevel > 2
    f(6) = figure('Position', [1 1 scrsz(3)*3/4 scrsz(4)/2]);
    subplot(1,2,1)
    barwitherr ([accuracy2.T1_Overall_stes(:,1) accuracy2.T1_Overall_stes(:,2)],[1 2 3 4],[accuracy2.T1_Overall_means(:,1) accuracy2.T1_Overall_means(:,2)])
    legend('valid','invalid')
    ylim([0 1])
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1]
    colormap(barmap);
    title('T1 overall')
    set(gca, 'XTick',[1 2 3 4],'XTickLabel',{'pp','pa','ap','aa' });
    
    subplot(1,2,2)
    barwitherr ([accuracy2.T2_Overall_stes(:,1) accuracy2.T2_Overall_stes(:,2)],[1 2 3 4],[accuracy2.T2_Overall_means(:,1) accuracy2.T2_Overall_means(:,2)])
    legend('valid','invalid')
    ylim([0 1])
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1];
    colormap(barmap);
    title('T2 overall')
    set(gca, 'XTick',[1 2 3 4],'XTickLabel',{'pp','pa','ap','aa' });
end

%% target axis: vertical, horizontal
accuracy3 = []; % rows: target axis [0,90] ; columns: '1-1','2-1','2-2','1-2'; pages: runs 1-10
targetAxis_Indx = [0,90];

for n = 1:length(df)
    runIndx = run == n;
    for k = 1:length(cueBlockOrder_Indx) % cueBlockOrder_Indx = [2 4 5 3]; cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'};
        dd = runIndx & cueBlockOrder == cueBlockOrder_Indx(k);
        if cueBlockOrder_Indx(k) == 2 || cueBlockOrder_Indx(k) == 4;
            interval = 1; TargetAxis = targetAxis(:,1);
            targetType = targetTypeT1;
        else interval = 2; TargetAxis = targetAxis(:,2);
            targetType = targetTypeT2;
        end
        for i = 1:2; % targetAxis = [0,90];
            condition_Indx =  dd & TargetAxis == targetAxis_Indx(i) ;
            % total correct / total present
            accuracy3.Discrim_all(i,k,n) = sum (targetType(condition_Indx) == response_all(condition_Indx))/ ...
                sum(ismember(targetType(condition_Indx),[1,2]) );
            % total correct / total correctly detected
            accuracy3.Discrim1_all(i,k,n) = sum (targetType(condition_Indx) == response_all(condition_Indx))/ ...
                sum (DetectTargetType(condition_Indx,interval) == DiscrimResponse_all(condition_Indx) );
            accuracy3.Detect_all(i,k,n) = sum (DetectTargetType(condition_Indx,interval) == DetectResponse_all(condition_Indx) ) /...
                numel(DetectTargetType(condition_Indx,interval));
            accuracy3.Overall_all(i,k,n) = sum ( response_correct(condition_Indx) == 1 ) / ...
                numel( response_correct(condition_Indx) );
            
        end
    end
end


%% calculate means and stes (vertical horizontal)

accuracy3.Discrim_means = nanmean(accuracy3.Discrim_all,3); % rows: [0 90]; columns: '1-1','2-1','2-2','1-2'
accuracy3.Discrim1_means = nanmean(accuracy3.Discrim1_all,3);
accuracy3.Detect_means = nanmean(accuracy3.Detect_all,3);
accuracy3.Overall_means = nanmean(accuracy3.Overall_all,3);

accuracy3.Discrim_stes = nanstd(accuracy3.Discrim_all,0,3)./sqrt(length(df));
accuracy3.Discrim1_stes = nanstd(accuracy3.Discrim1_all,0,3)./sqrt(length(df));
accuracy3.Detect_stes = nanstd(accuracy3.Detect_all,0,3)./sqrt(length(df));
accuracy3.Overall_stes = nanstd(accuracy3.Overall_all,0,3) ./ sqrt(length(df));

%% plot discrimination for vertical and horizontal
if plotLevel > 2
    f(7) = figure('Position', [1 1 scrsz(3)*3/4 scrsz(4)/2]);
    subplot(1,2,1)
    barwitherr ([accuracy3.Discrim1_stes(:,1) accuracy3.Discrim1_stes(:,2)],[1 2],[accuracy3.Discrim1_means(:,1) accuracy3.Discrim1_means(:,2)])
    legend('valid','invalid')
    ylim([0 1])
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1]
    colormap(barmap);
    title('T1 Discrim (total correct/total detected)')
    set(gca, 'XTick',[1 2],'XTickLabel',{'0','90' });
    
    subplot(1,2,2)
    barwitherr ([accuracy3.Discrim1_stes(:,3) accuracy3.Discrim1_stes(:,4)],[1 2],[accuracy3.Discrim1_means(:,3) accuracy3.Discrim1_means(:,4)])
    legend('valid','invalid')
    ylim([0 1])
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1];
    colormap(barmap);
    title('T2 Discrim (total correct/total detected)')
    set(gca, 'XTick',[1 2],'XTickLabel',{'0','90' });
end

%% plot detection for vertical and horizontal
if plotLevel > 2
    f(8) =  figure('Position', [1 1 scrsz(3)*3/4 scrsz(4)/2]);
    subplot(1,2,1)
    barwitherr ([accuracy3.Detect_stes(:,1) accuracy3.Detect_stes(:,2)],[1 2],[accuracy3.Detect_means(:,1) accuracy3.Detect_means(:,2)])
    legend('valid','invalid')
    ylim([0 1])
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1]
    colormap(barmap);
    title('T1 Detection')
    set(gca, 'XTick',[1 2],'XTickLabel',{'0','90' });
    
    subplot(1,2,2)
    barwitherr ([accuracy3.Detect_stes(:,3) accuracy3.Detect_stes(:,4)],[1 2],[accuracy3.Detect_means(:,3) accuracy3.Detect_means(:,4)])
    legend('valid','invalid')
    ylim([0 1])
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1];
    colormap(barmap);
    title('T2 Detection')
    set(gca, 'XTick',[1 2],'XTickLabel',{'0','90' });
end
%% plot overall accuracy for vertical and horizontal
if plotLevel > 2
    f(9) =  figure('Position', [1 1 scrsz(3)*3/4 scrsz(4)/2]);
    subplot(1,2,1)
    barwitherr ([accuracy3.Overall_stes(:,1) accuracy3.Overall_stes(:,2)],[1 2],[accuracy3.Overall_means(:,1) accuracy3.Overall_means(:,2)])
    legend('valid','invalid')
    ylim([0 1])
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1]
    colormap(barmap);
    title('T1 Overall')
    set(gca, 'XTick',[1 2],'XTickLabel',{'0','90' });
    
    subplot(1,2,2)
    barwitherr ([accuracy3.Overall_stes(:,3) accuracy3.Overall_stes(:,4)],[1 2],[accuracy3.Overall_means(:,3) accuracy3.Overall_means(:,4)])
    legend('valid','invalid')
    ylim([0 1])
    barmap=[0.7 0.7 0.7; 0.05 .45 0.1];
    colormap(barmap);
    title('T2 Overall')
    set(gca, 'XTick',[1 2],'XTickLabel',{'0','90' });
end

%% 2 targets: 2 axes x 2 orientations
accuracy4 = [];    % rows: runs = 1:10; columns: axis = [0 90]; pages: cueBlockOrder_Indx = [2 4 5 3]; 4thD: ori = [1 2];
sameOri = targetTypeT1 == targetTypeT2;
sameAxis = targetAxis (:,1) == targetAxis(:,2);
ori = [sameOri,~sameOri];
axis = [sameAxis,~sameAxis];

for n = 1:length(df)
    for o = 1:length(cueBlockOrder_Indx) % cueBlockOrder_Indx = [2 4 5 3]; cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'};
        if cueBlockOrder_Indx(o) == 2 || cueBlockOrder_Indx(o) == 4;
            interval = 1;
            targetType = targetTypeT1;
        else interval = 2;
            targetType = targetTypeT2;
        end
        for m = 1:2 % [sameOri, diffOri]
            for i = 1:2 % [sameAxis,diffAxis]
                condition_Indx  = run == n & cueBlockOrder == cueBlockOrder_Indx(o) & targetCondition == 2 & ...
                    ori(:,m) & axis(:,i);
                % total correct / total present
                accuracy4.Discrim_all(n,o,m,i) = sum (targetType(condition_Indx) == response_all(condition_Indx))/ ...
                    sum(ismember(targetType(condition_Indx),[1,2]) );
                % total correct / total correctly detected
                accuracy4.Discrim1_all(n,o,m,i) = sum (targetType(condition_Indx) == response_all(condition_Indx))/ ...
                    sum (DetectTargetType(condition_Indx,interval) == DiscrimResponse_all(condition_Indx) );
                accuracy4.Detect_all(n,o,m,i) = sum (DetectTargetType(condition_Indx,interval) == DetectResponse_all(condition_Indx) ) /...
                    numel(DetectTargetType(condition_Indx,interval));
                accuracy4.Overall_all(n,o,m,i) = sum ( response_correct(condition_Indx) == 1 ) / ...
                    numel( response_correct(condition_Indx) );
            end
        end
    end
end



% rows: cueBlockOrder = ['1-1','2-1','2-2','1-2']
% columns: orientation {same diff}
% pages: axis {same diff}
accuracy4.Discrim_means = squeeze (nanmean(accuracy4.Discrim_all,1));
accuracy4.Discrim1_means = squeeze (nanmean(accuracy4.Discrim1_all,1));
accuracy4.Detect_means = squeeze (nanmean(accuracy4.Detect_all,1));
accuracy4.Overall_means = squeeze (nanmean(accuracy4.Overall_all,1));

accuracy4.Discrim_stes = squeeze (nanstd(accuracy4.Discrim_all,0,1)./sqrt(length(df)));
accuracy4.Discrim1_stes = squeeze (nanstd(accuracy4.Discrim1_all,0,1)./sqrt(length(df)));
accuracy4.Detect_stes = squeeze (nanstd(accuracy4.Detect_all,0,1)./sqrt(length(df)));
accuracy4.Overall_stes = squeeze (nanstd(accuracy4.Overall_all,0,1) ./ sqrt(length(df)));

%% plot: 2 axes x 2 orientation
names = {'SOSA','SODA','DOSA','DODA'};

if plotLevel > 2
    for i = 1:2 % orientation {same diff}
        for k = 1:2 % axis {same diff}
            if i==1
                f(9+k)= figure;
            else
                f(9+k+2)= figure;
            end
            
            subplot(1,3,1)
            y = errorbar([ (accuracy4.Discrim1_means (1:2,i,k)') ;...
                (accuracy4.Discrim1_means (3:4,i,k)')],...
                [(accuracy4.Discrim1_stes(1:2,i,k)');...
                (accuracy4.Discrim1_stes(3:4,i,k)')],'.');
            ylim([0 1])
            set(y(2),'Color','r')
            set(gca,'XTick',[1 2])
            set(gca,'XTickLabel',{'T1','T2'});
            ylabel('accuracy')
            legend('valid','invalid','Location','SouthEast');
            if i == 1
                title([names{k},' - discrimination'])
            else
                title ([names{k+2},' - discrimination'])
            end
            
            subplot (1,3,2)
            y = errorbar([ (accuracy4.Detect_means (1:2,i,k)') ;...
                (accuracy4.Detect_means (3:4,i,k)')],...
                [(accuracy4.Detect_stes(1:2,i,k)');...
                (accuracy4.Detect_stes(3:4,i,k)')],'.');
            ylim([0 1])
            set(y(2),'Color','r')
            set(gca,'XTick',[1 2])
            set(gca,'XTickLabel',{'T1','T2'});
            ylabel('accuracy')
            legend('valid','invalid','Location','SouthEast');
            if i == 1
                title([names{k},' - detection'])
            else
                title ([names{k+2},' - detection'])
            end
            
            subplot(1,3,3)
            y = errorbar([ (accuracy4.Overall_means (1:2,i,k)') ;...
                (accuracy4.Overall_means (3:4,i,k)')],...
                [(accuracy4.Overall_stes(1:2,i,k)');...
                (accuracy4.Overall_stes(3:4,i,k)')],'.');
            ylim([0 1])
            set(y(2),'Color','r')
            set(gca,'XTick',[1 2])
            set(gca,'XTickLabel',{'T1','T2'});
            ylabel('accuracy')
            legend('valid','invalid','Location','SouthEast');
            if i == 1
                title([names{k},' - overall'])
            else
                title ([names{k+2},' - overall'])
            end
            
            turnallwhite
            
        end
    end
end

%% turn figs white
turnallwhite

%% Save analysis files and figures

% save analysis
if saveFile
    save(sprintf('%s/%s.mat', analysisDir, analysisFile), 'accuracy',...
        'accuracy2','accuracy3','accuracy4','responseData_all',...
        'responseData_labels')
end

% save figs
if saveFigs
    if plotLevel == 1
        rd_saveAllFigs(f, {'minAll'},analysisFile,figDir);
    elseif plotLevel == 2
        rd_saveAllFigs([], {'minAll','Discrim2','Detect2'},analysisFile,figDir);
    elseif plotLevel == 3
        rd_saveAllFigs(f, {'all','Detect1','Discrim1','Discrim2',...
            'Detect2','overall','DiscrimAxis','DetectAxis','OverallAxis','SOSA','SODA','DOSA','DODA'},analysisFile,figDir);
    end
end

%% Store io info
io.analysisFile = analysisFile;
io.figDir = figDir;

%%
% figure
% y = errorbar([Pc.Discrim_all_mean;Pc.Detect_all_mean],[Pc.Discrim_all_ste;Pc.Detect_all_ste],'.');
% ylim([0 1])
% set(gca,'XTickLabel',{'','discrimination','','','','','detection',''});
% legend(y,{'T1','T2'},'Location','NorthWest')
% ylabel('Accuracy')
% title('overall accuracy')

end


