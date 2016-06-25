% rd_analyzeTADetectDiscrimGroupBehav.m

%% setup
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';

subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
    'R0861_20150813', 'R0504_20150805', 'R0983_20150813', ...
    'R0898_20150828', 'R0436_20150904', 'R1018_20151118', ...
    'R1019_20151118','R1021_20151120','R1026_20151211', ...
    'R0852_20151211','R1027_20151216','R1028_20151216',...
    'R1029_20151222'}; % N=16

nSubjects = numel(subjects);

%% load data
for iSubject = 1:nSubjects
    sessionDir = subjects{iSubject}; 
    behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);
    behavFile = dir(sprintf('%s/*.mat', behavDir));
    b = load(sprintf('%s/%s', behavDir, behavFile.name));
    behav(iSubject) = behavior(b); % update behav with more info
end

%% organize data from Sirui (accuracy field)
% accuracy rows: '1-1','2-1','2-2','1-2'
measures = {'Discrim_means','Discrim1_means','Detect_means','Overall_means',...
    'Hit_means','FA_means','Miss_means','CR_means','dprime_means'};
nM = numel(measures);

for iSubject = 1:nSubjects
    for iM = 1:nM
        m = measures{iM};
        groupDataS.(m)(:,iSubject) = behav(iSubject).accuracy.(m);
    end
end

for iM = 1:nM
    m = measures{iM};
    groupMeanS(:,iM) = mean(groupDataS.(m),2);
    groupSteS(:,iM) = std(groupDataS.(m),0,2)/sqrt(nSubjects);
end

%% plot
figure('color','w')
for iM = 1:nM
    subplot(1,nM,iM)
    m = measures{iM};
    p1 = errorbar(reshape(groupMeanS(:,iM),[2,2])', reshape(groupSteS(:,iM),[2,2])');
    set(p1(2),'color','r')
    xlim([.5 2.5])
    ylim([0 1])
    title(m)
    set(gca,'XTick',[1 2])
    box off
end
ylim([0 3])
legend('valid','invalid')

%% reanalyze data
targets = unique(behav(1).responseTarget(behav(1).responseTarget~=0));
cueValidities = unique(behav(1).cueValidity(behav(1).cueValidity~=0));
cueValidities = sort(cueValidities,'descend'); % 1=valid, -1=invalid

groupDataAll = [];
for iS = 1:nSubjects
    for iT = 1:numel(targets)
        wT = behav(iS).responseTarget==targets(iT);
        for iV = 1:numel(cueValidities)
            wV = behav(iS).cueValidity==cueValidities(iV);
            w = wT & wV;
            groupDataAll(iS).detectHMFC{iV,iT} = behav(iS).detectHMFC(w,:);
            groupDataAll(iS).discrimCI{iV,iT} = behav(iS).discrimCI(w,:);
            groupDataAll(iS).acc{iV,iT} = behav(iS).acc(w,:);
        end
    end
end

for iS = 1:nSubjects
    detect = groupDataAll(iS).detectHMFC;
    discrim = groupDataAll(iS).discrimCI;
    acc = groupDataAll(iS).acc;
    for iT = 1:numel(targets)
        for iV = 1:numel(cueValidities)
            presentResponse = any(discrim{iV,iT},2);
            groupData.discrim(iV,iT,iS) = nanmean(discrim{iV,iT}(:,1));
            groupData.discrim1(iV,iT,iS) = nanmean(discrim{iV,iT}(presentResponse,1));
            groupData.hit(iV,iT,iS) = nanmean(detect{iV,iT}(:,1));
            groupData.miss(iV,iT,iS) = nanmean(detect{iV,iT}(:,2));
            groupData.fa(iV,iT,iS) = nanmean(detect{iV,iT}(:,3));
            groupData.cr(iV,iT,iS) = nanmean(detect{iV,iT}(:,4));
            groupData.overallAcc(iV,iT,iS) = nanmean(acc{iV,iT});
        end
    end
end

% calculate dprime
h = groupData.hit;
fa = groupData.fa;
h(h==1) = .99;
fa(fa==0) = .01;
groupData.dprimeDetect = norminv(h) - norminv(fa);
groupData.critDetect = -.5 * (norminv(h) + norminv(fa));

%% group summary
measures = fields(groupData);
nM = numel(measures);
for iM = 1:nM
    m = measures{iM};
    groupMean.(m) = mean(groupData.(m),3);
    groupSte.(m) = std(groupData.(m),0,3)./sqrt(nSubjects);
end

%% plot
figure('color','w')
for iM = 1:nM
    subplot(1,nM,iM)
    m = measures{iM};
    p1 = errorbar(groupMean.(m)', groupSte.(m)');
    set(p1(2),'color','r')
    xlim([.5 2.5])
    switch m
        case 'dprimeDetect'
            ylim([0 3])
        case 'critDetect'
            ylim([-1 1])
            hold on
            plot([.5 2.5],[0 0],'--k')
        otherwise
            ylim([0 1])
    end
    title(m)
    set(gca,'XTick',[1 2])
    box off
end
legend('valid','invalid')

indivM = {'overallAcc','discrim1','dprimeDetect','critDetect'};
figure('Color','w','Position',[0 0 1200 800])
for iS = 1:nSubjects
    for iM = 1:numel(indivM)
        row = ceil(iS/2);
        if mod(iS,2)
            col = iM;
        else
            col = iM+5;
        end
        subplot(8,9,(row-1)*9+col)
        m = indivM{iM};
        p1 = plot(groupData.(m)(:,:,iS)');
        set(p1(2),'color','r')
        xlim([.5 2.5])
        switch m
            case 'dprimeDetect'
                ylim([0 3])
            case 'critDetect'
                ylim([-1 1])
                hold on
                plot([.5 2.5],[0 0],'--k')
            otherwise
                ylim([0 1])
        end
        if row==1
            title(m)
        else
%             set(gca,'XTickLabel','')
%             set(gca,'YTickLabel','')
        end
        if col==1 || col==6
%             ylabel(sprintf('S%d',iS))
        end
        set(gca,'XTick',[1 2])
        box off
    end
end
legend('valid','invalid')

%% stats
% valid vs. invalid
for iM = 1:numel(indivM)
    m = indivM{iM};
    vi.(m) = squeeze((groupData.(m)(1,:,:)-groupData.(m)(2,:,:)))';
    [h, pstat.(m), ci, stat.(m)] = ttest(vi.(m));
end




