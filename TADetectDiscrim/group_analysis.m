function [accuracy,accuracy2,responseData_group] = group_analysis (subjects,runs,dates, plotLevel)
rootDir = pathToExpt;
nSubj = length(subjects);
responseData_group = [];

if ~exist('dates','var')
    date = [];
end

if ~exist('plotLevel','var')
    plotLevel = 3; % 1 plot everything --> 3 plot minimal
end


for n = 1:nSubj;
    subject = subjects{n};
    date = dates{n};
    dir = sprintf('%s/analysis/%s', rootDir, subject);
    if numel(runs)>1
        analysisFile = sprintf('%s%s_taDetectDiscrim_%s%s%s%s',subject,date,num2str(runs(1)),'_',num2str(runs(end)));
    else
        analysisFile = sprintf('%s_%s_taDetectDiscrim_%s%s',subject,date,num2str(runs)');
    end
    analysis = load (sprintf('%s/%s.mat', dir, analysisFile));
    temp = analysis.responseData_all;
    responseData_group = [responseData_group;temp];
    
end

%% extract block order from responseData_group
run = responseData_group(:,1);
cueBlockOrder = responseData_group(:,4);
targetCondition = responseData_group(:,6);
targetTypeT1 = responseData_group(:,7);
targetTypeT2 = responseData_group(:,8);
response_all = responseData_group(:,10);
response_correct = responseData_group(:,11);
targetAxis = responseData_group(:,13:14);

% convert target type and response data for computing detection rate
DetectTargetType = responseData_group(:,7:8);
DetectTargetType (DetectTargetType == 2) = 1;
DetectResponse_all = response_all;
DetectResponse_all (DetectResponse_all  == 2) = 1;
DetectResponse_all (DetectResponse_all  == 3) = 0;

DiscrimResponse_all = response_all;
DiscrimResponse_all (DiscrimResponse_all  == 2) = 1;

OverallResponse_all = response_all;
OverallResponse_all ( OverallResponse_all == 3) = 0;



%% For each run calculate detection and discrimination accuracy for:
% postcue T1 valid & invalid; postcue T2 valid & invalid; overall accuracy

% blockNames = {'blank','fast-left'}; % fast-left
% attBlockNames = {'no-att','att-right'}; % att-right
% targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres','abs-abs'};
% cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'}; % 2-1 = cueT2,postcueT1

cueBlockOrder_Indx = [2 4 5 3];
accuracy =[]; % rows: '1-1','2-1','2-2','1-2'; columns: runs 1-10
accuracy.missedTrials = zeros(1,length(runs)); % columns: runs 1-10

for n = 1:length(runs)
    runIndx = run == n;
    for k = 1:length(cueBlockOrder_Indx) % cueBlockOrder_Indx = [2 4 5 3]; cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'};
        condition_Indx = runIndx & cueBlockOrder == cueBlockOrder_Indx(k);
        if cueBlockOrder_Indx(k) == 2 || cueBlockOrder_Indx(k) == 4;
            interval = 1;
            targetType = targetTypeT1;
        else interval = 2;
            targetType = targetTypeT2;
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
accuracy.Discrim_stes = nanstd(accuracy.Discrim_all,0,2)./sqrt(length(runs));
accuracy.Discrim1_stes = nanstd(accuracy.Discrim1_all,0,2)./sqrt(length(runs));
accuracy.Detect_stes = nanstd(accuracy.Detect_all,0,2)./sqrt(length(runs));
accuracy.Overall_stes = nanstd(accuracy.Overall_all,0,2) ./ sqrt(length(runs));
accuracy.Hit_stes = nanstd(accuracy.Hit_all,0,2) ./ sqrt(length(runs));
accuracy.FA_stes = nanstd(accuracy.FA_all,0,2)./ sqrt(length(runs));
accuracy.Miss_stes = nanstd(accuracy.Miss_all,0,2)./ sqrt(length(runs));
accuracy.CR_stes = nanstd(accuracy.CR_all,0,2)./ sqrt(length(runs));
accuracy.dprime_stes = nanstd(accuracy.dprime,0,2) ./ sqrt(length(runs));
all_stes = [accuracy.Hit_stes,accuracy.FA_stes,accuracy.Miss_stes,accuracy.CR_stes];



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
        legend('valid','invalid','Location','SouthEast');
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
        f(1) = figure('Position', [1 scrsz(4) scrsz(3)*(3/5) scrsz(4)/2]);
        
        subplot(1,3,1)
        y = errorbar([ (accuracy.Detect_means(1:2)');(accuracy.Detect_means(3:4)')],[(accuracy.Detect_stes(1:2)');(accuracy.Detect_stes(3:4)')],'.');
        ylim([0 1])
        % set(y, 'MarkerSize', 20)
        set(y(2),'Color','r')
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'});
        % set(gca,'XTickLabel',{'','T1','','','','','T2',''});
        ylabel('accuracy')
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
        error('plotLevel not recognized')
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

%% pp pa ap aa
accuracy2 = []; % rows: pp,pa,ap,aa; columns: '1-1','2-1','2-2','1-2'; pages: runs 1-10

for n = 1:length(runs)
    runIndx = run == n;
    for k = 1:length(cueBlockOrder_Indx) % cueBlockOrder_Indx = [2 4 5 3]; cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'};
        dd = runIndx & cueBlockOrder == cueBlockOrder_Indx(k);
        for i = 2:5; % targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres','abs-abs'};
            condition_Indx =  dd & targetCondition == i ;
            if cueBlockOrder_Indx(k) == 2 || cueBlockOrder_Indx(k) == 4;
                interval = 1;
                targetType = targetTypeT1;
            else interval = 2;
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

accuracy2.Discrim_stes = nanstd(accuracy2.Discrim_all,0,3)./sqrt(length(runs));
accuracy2.Discrim1_stes = nanstd(accuracy2.Discrim1_all,0,3)./sqrt(length(runs));
accuracy2.Detect_stes = nanstd(accuracy2.Detect_all,0,3)./sqrt(length(runs));
accuracy2.Overall_stes = nanstd(accuracy2.Overall_all,0,3) ./ sqrt(length(runs));

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
if plotLevel==1
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

if plotLevel < 3
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

%% plot detection for pp ap ap aa
if plotLevel < 3
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
if plotLevel==1
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

for n = 1:length(runs)
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

accuracy3.Discrim_stes = nanstd(accuracy3.Discrim_all,0,3)./sqrt(length(runs));
accuracy3.Discrim1_stes = nanstd(accuracy3.Discrim1_all,0,3)./sqrt(length(runs));
accuracy3.Detect_stes = nanstd(accuracy3.Detect_all,0,3)./sqrt(length(runs));
accuracy3.Overall_stes = nanstd(accuracy3.Overall_all,0,3) ./ sqrt(length(runs));

%% plot discrimination for vertical and horizontal
if plotLevel == 1
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
if plotLevel == 1
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
if plotLevel == 1
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

for n = 1:length(runs)
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

accuracy4.Discrim_stes = squeeze (nanstd(accuracy4.Discrim_all,0,1)./sqrt(length(runs)));
accuracy4.Discrim1_stes = squeeze (nanstd(accuracy4.Discrim1_all,0,1)./sqrt(length(runs)));
accuracy4.Detect_stes = squeeze (nanstd(accuracy4.Detect_all,0,1)./sqrt(length(runs)));
accuracy4.Overall_stes = squeeze (nanstd(accuracy4.Overall_all,0,1) ./ sqrt(length(runs)));

%% plot: 2 axes x 2 orientation
names = {'SOSA','SODA','DOSA','DODA'};

if plotLevel == 1
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

% saveFile = 1;
% saveFigs = 1;

% % save analysis
% if saveFile
%     save(sprintf('%s/%s.mat', dataDir, analysisFile), 'accuracy',...
%         'accuracy2','responseData_group','responseData_labels')
% end
%
% % save figs
% if saveFigs
%     if plotLevel == 1;
%     rd_saveAllFigs(f, {'all','Detect1','Discrim1','Discrim2',...
%         'Detect2','overall','DiscrimAxis','DetectAxis','OverallAxis'},...
%         analysisFile,figDir);
%     elseif plotLevel == 2;
%     rd_saveAllFigs([], {'minAll','Discrim2','Detect2'},analysisFile,figDir);
%     elseif plotLevel == 3;
%     rd_saveAllFigs(f, {'minAll'},analysisFile,figDir);
%     end
% end


%%
% figure
% y = errorbar([Pc.Discrim_all_mean;Pc.Detect_all_mean],[Pc.Discrim_all_ste;Pc.Detect_all_ste],'.');
% ylim([0 1])
% set(gca,'XTickLabel',{'','discrimination','','','','','detection',''});
% legend(y,{'T1','T2'},'Location','NorthWest')
% ylabel('Accuracy')
% title('overall accuracy')

end



