function [accuracy, stim] = rd_analyzeTADetectDiscrimOneRun(dataDir, stimDir, stimfile, run, plotLevel)

dataFile = dir(sprintf('%s/*%s%d.mat', dataDir, stimfile, run));
dd = load(sprintf('%s/%s', dataDir, dataFile(end).name));
stim = load(sprintf('%s/%s%d.mat', stimDir, stimfile, run));

trialCount = length(stim.p.blockOrder);  
respSecs = stim.p.respDur;
feedbackDur = stim.p.feedbackDur;
refreshRate = stim.p.refrate;  %(frames)
blockLength = stim.p.blockDur*refreshRate; %(frames)
% respTime = blockLength-(respSecs+feedbackDur)*refreshRate; %frames to respond period
respTime = (stim.p.targetLeadTime + stim.p.targetSOA + stim.p.cueTargetSOA)*refreshRate;
keyCodes = stim.p.keyCodes;
order = stim.order;
responseOption = stim.p.responseOption;

if isfield(stim.stimulus.target, 'catchTrials')
    catchTrials = stim.stimulus.target.catchTrials;
else
    catchTrials = false;
end

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
n = 1;
df = 1:n; % placeholder for list of data files

[responseData_all, responseData_labels] = ...
    sl_responseDiscrimData(respTime,trialCount,...
    respSecs,refreshRate, blockLength, keyCodes, dd.response, ...
    order, n, itiSeq, jitSeq, stim.p.jitter);

%% extract block order from responseData_all
run = responseData_all(:,1);
cueBlockOrder = responseData_all(:,4);
targetCondition = responseData_all(:,6);
targetTypeT1 = responseData_all(:,7);
targetTypeT2 = responseData_all(:,8);
response_all = responseData_all(:,10);
response_correct = responseData_all(:,11);
targetAxis = responseData_all(:,13:14);
targetPedestal = responseData_all(:,15:16);

if isfield(stim.order,'posBlockOrder')
    targetPos = stim.order.posBlockOrder';
    targetPosType = targetPos;
    targetPosType(targetPos>=1 & targetPos<=4) = 1;
    targetPosType(targetPos>=5 & targetPos<=8) = 2;
end

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
cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'}; % 2-1 = cueT2,postcueT1

cueBlockOrder_Indx = [2 4 5 3]; 
accuracy =[]; % rows: '1-1','2-1','2-2','1-2'; columns: runs 1-10
accuracy.missedTrials = zeros(1,length(df)); % columns: runs 1-10


runIndx = 1;
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
                targetType = targetPosType(:,1);
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
                targetType = targetPosType(:,2);
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
    
    accuracy.missedTrials(n) = accuracy.missedTrials(n) + sum(isnan(response_all(condition_Indx)) );
end

%% Calculate target type 1 valid & invalid; target type 2 valid and invalid
for k = 1:numel(cueBlockOrder_Indx) %'1-1','2-1','2-2','1-2'
    cueType = cueBlockOrder_Indx(k);
    cueName = cueBlockNames{cueType};
    target = str2num(cueName(end));
    switch target
        case 1
            targetTypeOrder = targetTypeT1;
            pedestalOrder = targetPedestal(:,1);
        case 2
            targetTypeOrder = targetTypeT2;
            pedestalOrder = targetPedestal(:,2);
    end
    for tt = 1:2 % target type, e.g. decrement vs. increment
        w = cueBlockOrder==cueType & targetTypeOrder==tt;
        temp_correct = response_correct(w); 
        temp_correct(temp_correct==-1) = 0;
        ttacc{k,tt} = temp_correct;
        ttacc_mean(k,tt) = mean(temp_correct);
    end
    if catchTrials
        tt = 0; % target absent catch trials
        w = cueBlockOrder==cueType & targetTypeOrder==tt;
        ttcatch{k,1} = OverallResponse_all(w);
    end
    for tped = 1:2 % pedestal, e.g. low contrast vs. high contrast
        w = cueBlockOrder==cueType & pedestalOrder==tped;
        temp_correct = response_correct(w); 
        temp_correct(temp_correct==-1) = 0;
        tpedacc{k,tped} = temp_correct;
        tpedacc_mean(k,tped) = mean(temp_correct);
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
all_stes = [accuracy.Hit_stes,accuracy.FA_stes,accuracy.Miss_stes,accuracy.CR_stes];

% for single run
accuracy.targetTypeAccAll = ttacc; 
accuracy.targetTypeAcc = ttacc_mean; 
if catchTrials
    accuracy.catchTrialRespAll = ttcatch;
end
accuracy.targetPedestalAccAll = tpedacc;
accuracy.targetPedestalAcc = tpedacc_mean;

%% plot
switch plotLevel
    case 1
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
        
    case {2,3}
        scrsz=get(0,'ScreenSize');
%         f(1) = figure('Position', [1 scrsz(4) scrsz(3)*(3/5) scrsz(4)/2]);
        f(1) = figure('Position', [1 600 800 400]);
        
        subplot(1,3,1)
        y = errorbar([ (accuracy.Detect_means(1:2)');(accuracy.Detect_means(3:4)')],[(accuracy.Detect_stes(1:2)');(accuracy.Detect_stes(3:4)')],'.');
        ylim([0 1])
        % set(y, 'MarkerSize', 20)
        set(y(2),'Color','r')
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'});
        % set(gca,'XTickLabel',{'','T1','','','','','T2',''});
        ylabel('accuracy')
        legend('valid','invalid','Location','SouthEast');
        title('detection')
        
        subplot(1,3,2)
        y = errorbar([ (accuracy.Discrim1_means(1:2)'); (accuracy.Discrim1_means(3:4)')],[(accuracy.Discrim1_stes(1:2)');(accuracy.Discrim1_stes(3:4)')],'.');
        ylim([0 1])
        set(y(2),'Color','r')
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'});
        ylabel('accuracy')
        title('discrimination(total correct/total detected)')
        
        subplot(1,3,3)
        y = errorbar([ (accuracy.Overall_means(1:2)');(accuracy.Overall_means(3:4)')],[(accuracy.Overall_stes(1:2)');(accuracy.Overall_stes(3:4)')],'.');
        ylim([0 1])
        set(y(2),'Color','r')
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'});
        ylabel('accuracy')
        title('overall')
    otherwise
        % no plots
end

%% plot (hit, fa, miss, cr)
if plotLevel==1
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

% %% pp pa ap aa
% accuracy2 = []; % rows: pp,pa,ap,aa; columns: '1-1','2-1','2-2','1-2'; pages: runs 1-10
% 
% for n = 1:length(df) 
%     runIndx = run == n;
%     for k = 1:length(cueBlockOrder_Indx) % cueBlockOrder_Indx = [2 4 5 3]; cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'}; 
%         dd = runIndx & cueBlockOrder == cueBlockOrder_Indx(k);
%         for i = 2:5; % targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres','abs-abs'};
%             condition_Indx =  dd & targetCondition == i ;
%             if cueBlockOrder_Indx(k) == 2 || cueBlockOrder_Indx(k) == 4;
%                  interval = 1;
%                  targetType = targetTypeT1;
%             else interval = 2;
%                  targetType = targetTypeT2;
%             end
%               % total correct / total present
%              accuracy2.Discrim_all(i-1,k,n) = sum (targetType(condition_Indx) == response_all(condition_Indx))/ ...
%                  sum(ismember(targetType(condition_Indx),[1,2]) );
%              % total correct / total correctly detected
%              accuracy2.Discrim1_all(i-1,k,n) = sum (targetType(condition_Indx) == response_all(condition_Indx))/ ...
%                  sum (DetectTargetType(condition_Indx,interval) == DiscrimResponse_all(condition_Indx) );
%              accuracy2.Detect_all(i-1,k,n) = sum (DetectTargetType(condition_Indx,interval) == DetectResponse_all(condition_Indx) ) /...
%                  numel(DetectTargetType(condition_Indx,interval));
%              accuracy2.Overall_all(i-1,k,n) = sum ( response_correct(condition_Indx) == 1 ) / ...
%                  numel( response_correct(condition_Indx) );
%              
%          end
%      end
% end
%  
% %% calculate means and stes (pp pa ap aa)
% accuracy2.Discrim_means = nanmean(accuracy2.Discrim_all,3); % rows:  pp,pa,ap,aa; columns: target type
% accuracy2.Discrim1_means = nanmean(accuracy2.Discrim1_all,3);
% accuracy2.Detect_means = nanmean(accuracy2.Detect_all,3);
% accuracy2.Overall_means = nanmean(accuracy2.Overall_all,3);
% 
% accuracy2.Discrim_stes = nanstd(accuracy2.Discrim_all,0,3)./sqrt(length(df));
% accuracy2.Discrim1_stes = nanstd(accuracy2.Discrim1_all,0,3)./sqrt(length(df));
% accuracy2.Detect_stes = nanstd(accuracy2.Detect_all,0,3)./sqrt(length(df));
% accuracy2.Overall_stes = nanstd(accuracy2.Overall_all,0,3) ./ sqrt(length(df));
% 
% accuracy2.T1_Discrim_means = accuracy2.Discrim_means(1:2,1:2); % rows: pp,pa; columns: T1 valid, invalid; 
% accuracy2.T1_Discrim_stes = accuracy2.Discrim_stes(1:2,1:2);
% accuracy2.T2_Discrim_means = [accuracy2.Discrim_means(1,3:4) ; accuracy2.Discrim_means(3,3:4)]; % rows: pp,ap columns: T2 valid, invalid;
% accuracy2.T2_Discrim_stes = [accuracy2.Discrim_stes(1,3:4) ; accuracy2.Discrim_stes(3,3:4)];
% 
% accuracy2.T1_Discrim1_means = accuracy2.Discrim1_means(1:2,1:2); % rows: pp,pa; columns:  T1 valid, invalid;
% accuracy2.T1_Discrim1_stes = accuracy2.Discrim1_stes(1:2,1:2);
% accuracy2.T2_Discrim1_means = [accuracy2.Discrim1_means(1,3:4) ; accuracy2.Discrim1_means(3,3:4)];% rows: pp,ap columns: T2 valid, invalid;
% accuracy2.T2_Discrim1_stes = [accuracy2.Discrim1_stes(1,3:4) ; accuracy2.Discrim1_stes(3,3:4)];
% 
% accuracy2.T1_Detect_means = accuracy2.Detect_means(:,1:2);
% accuracy2.T1_Detect_stes = accuracy2.Detect_stes(:,1:2);
% accuracy2.T2_Detect_means = accuracy2.Detect_means(:,3:4);
% accuracy2.T2_Detect_stes = accuracy2.Detect_stes(:,3:4);
% 
% accuracy2.T1_Overall_means = accuracy2.Overall_means(:,1:2);
% accuracy2.T1_Overall_stes = accuracy2.Overall_stes(:,1:2);
% accuracy2.T2_Overall_means = accuracy2.Overall_means(:,3:4);
% accuracy2.T2_Overall_stes = accuracy2.Overall_stes(:,3:4);
% 
% 
% %% plot discrimination pp pa/pp ap for T1 and T2
% if plotLevel==1
%     f(3) = figure('Position', [1 scrsz(4) scrsz(3)*3/4 scrsz(4)/2]);
%     subplot(1,2,1)
%     ylim([0 1])
%     barwitherr ([accuracy2.T1_Discrim_stes(:,1) accuracy2.T1_Discrim_stes(:,2)],[1 2],[accuracy2.T1_Discrim_means(:,1) accuracy2.T1_Discrim_means(:,2)])
%     legend('valid','invalid','Location','best')
%     barmap=[0.7 0.7 0.7; 0.05 .45 0.1]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1]
%     colormap(barmap);
%     title('T1 Discrim (total correct/total present)')
%     set(gca, 'XTick',[1 2],'XTickLabel',{'pp','pa' });
%     
%     subplot(1,2,2)
%     barwitherr ([accuracy2.T2_Discrim_stes(:,1) accuracy2.T2_Discrim_stes(:,2)],[1 2],[accuracy2.T2_Discrim_means(:,1) accuracy2.T2_Discrim_means(:,2)])
%     legend('valid','invalid','Location','best')
%     ylim([0 1])
%     barmap=[0.7 0.7 0.7; 0.05 .45 0.1];
%     colormap(barmap);
%     title('T2 Discrim (total correct/total present)')
%     set(gca, 'XTick',[1 2],'XTickLabel',{'pp','ap' });
% end
% 
% if plotLevel < 3
%     f(4) = figure('Position', [1 scrsz(1) scrsz(3)*3/4 scrsz(4)/2]);
%     subplot(1,2,1)
%     
%     barwitherr ([accuracy2.T1_Discrim1_stes(:,1) accuracy2.T1_Discrim1_stes(:,2)],[1 2],[accuracy2.T1_Discrim1_means(:,1) accuracy2.T1_Discrim1_means(:,2)])
%     legend('valid','invalid','Location','NorthEast')
%     ylim([0 1])
%     barmap=[0.7 0.7 0.7; 0.05 .45 0.1]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1]
%     colormap(barmap);
%     title('T1 Discrim (total correct/total detected)')
%     set(gca, 'XTick',[1 2],'XTickLabel',{'pp','pa' });
%     
%     subplot(1,2,2)
%     barwitherr ([accuracy2.T2_Discrim1_stes(:,1) accuracy2.T2_Discrim1_stes(:,2)],[1 2],[accuracy2.T2_Discrim1_means(:,1) accuracy2.T2_Discrim1_means(:,2)])
%     legend('valid','invalid','Location','NorthEast')
%     ylim([0 1])
%     barmap=[0.7 0.7 0.7; 0.05 .45 0.1];
%     colormap(barmap);
%     title('T2 Discrim (total correct/total detected)')
%     set(gca, 'XTick',[1 2],'XTickLabel',{'pp','ap' });
% end
% 
% %% plot detection for pp ap ap aa
% if plotLevel < 3
%     f(5) = figure('Position', [1 1 scrsz(3)*3/4 scrsz(4)/2]);
%     subplot(1,2,1)
%     barwitherr ([accuracy2.T1_Detect_stes(:,1) accuracy2.T1_Detect_stes(:,2)],[1 2 3 4],[accuracy2.T1_Detect_means(:,1) accuracy2.T1_Detect_means(:,2)])
%     legend('valid','invalid')
%     ylim([0 1])
%     barmap=[0.7 0.7 0.7; 0.05 .45 0.1]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1]
%     colormap(barmap);
%     title('T1 detection')
%     set(gca, 'XTick',[1 2 3 4],'XTickLabel',{'pp','pa','ap','aa' });
%     
%     subplot(1,2,2)
%     barwitherr ([accuracy2.T2_Detect_stes(:,1) accuracy2.T2_Detect_stes(:,2)],[1 2 3 4],[accuracy2.T2_Detect_means(:,1) accuracy2.T2_Detect_means(:,2)])
%     legend('valid','invalid')
%     ylim([0 1])
%     barmap=[0.7 0.7 0.7; 0.05 .45 0.1];
%     colormap(barmap);
%     title('T2 detection')
%     set(gca, 'XTick',[1 2 3 4],'XTickLabel',{'pp','pa','ap','aa' });
% end
% 
% %% plot overall accuracy for pp pa ap aa
% if plotLevel==1
%     f(6) = figure('Position', [1 1 scrsz(3)*3/4 scrsz(4)/2]);
%     subplot(1,2,1)
%     barwitherr ([accuracy2.T1_Overall_stes(:,1) accuracy2.T1_Overall_stes(:,2)],[1 2 3 4],[accuracy2.T1_Overall_means(:,1) accuracy2.T1_Overall_means(:,2)])
%     legend('valid','invalid')
%     ylim([0 1])
%     barmap=[0.7 0.7 0.7; 0.05 .45 0.1]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1]
%     colormap(barmap);
%     title('T1 overall')
%     set(gca, 'XTick',[1 2 3 4],'XTickLabel',{'pp','pa','ap','aa' });
%     
%     subplot(1,2,2)
%     barwitherr ([accuracy2.T2_Overall_stes(:,1) accuracy2.T2_Overall_stes(:,2)],[1 2 3 4],[accuracy2.T2_Overall_means(:,1) accuracy2.T2_Overall_means(:,2)])
%     legend('valid','invalid')
%     ylim([0 1])
%     barmap=[0.7 0.7 0.7; 0.05 .45 0.1];
%     colormap(barmap);
%     title('T2 overall')
%     set(gca, 'XTick',[1 2 3 4],'XTickLabel',{'pp','pa','ap','aa' });
% end

