% rd_analyzeTADetectDiscrimGroupBehav.m

%% setup
% exptDir = pathToExpt;
% exptDir = '/Local/Users/denison/Data/TANoise/Behavior';
% exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/Behavior';
% exptDir = '/Local/Users/denison/Data/TADetectDiscrim/Behavior';
% exptDir = '/Local/Users/denison/Data/TA2/Behavior';
% exptDir = pathToTA2('Behavior');
exptDir = pathToTANoise('Behavior');

% subjects = {'rd','lr','mj','af','xw'};
% startRuns = [211, 211, 221, 221, 221];

% subjects = {'jpnoise'};
% startRuns = 16;

subjects = {'R0817_20171212','R0817_20171213', ...
        'R0898_20180112','R0898_20180116',...
        'R0959_20180219','R0959_20180306',...
        'R0983_20180111','R0983_20180112', ...
        'R1021_20180208','R1021_20180212',...
        'R1103_20180213','R1103_20180215',...
        'R1187_20180105','R1187_20180108',...
        'R1452_20190717','R1452_20190718',...
        'R1507_20190702','R1507_20190705',...
        'R1535_20190717','R1535_20190718',...
        }; % TANoise
startRuns = repmat(1,numel(subjects));

% subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
%     'R0861_20150813', 'R0504_20150805', 'R0983_20150813', ...
%     'R0898_20150828', 'R0436_20150904', 'R1018_20151118', ...
%     'R1019_20151118','R1021_20151120','R1026_20151211', ...
%     'R0852_20151211','R1027_20151216','R1028_20151216',...
%             'R1029_20151222'}; % N=16 TADetectDiscrim
% % subjects = subjects([1 2 4 5 7 8 10 12 14 16]);
% startRuns = repmat(1,numel(subjects));

%     'R0890_20181121', ...
% subjects = {'R0817_20181120', 'R0817_20190625',...
%     'R0959_20181128', 'R0959_20190703',...
%     'R1103_20181121', 'R1103_20190710', ...
%     'R1187_20181119', 'R1187_20190703',...
%     'R1373_20181128', 'R1373_20190708',...
%     'R1452_20181119', 'R1452_20190711',...
%     'R1507_20190621', 'R1507_20190627',...
%     'R1535_20190708', 'R1535_20190711'}; % N=8 TA2
% subjects = subjects([1 2 4 5 7 8 10 12 14 16]);
startRuns = repmat(1,numel(subjects));

nSubjects = numel(subjects);

%% load data
for iSubject = 1:nSubjects
    sessionDir = subjects{iSubject}; 
%     behavDir = sprintf('%s/analysis/%s', exptDir, sessionDir);
    behavDir = sprintf('%s/%s/analysis', exptDir, sessionDir);
    behavFile = dir(sprintf('%s/*%d*.mat', behavDir, startRuns(iSubject)));
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
            groupDataAll(iS).rt{iV,iT} = behav(iS).rt(w,:);
        end
    end
end

for iS = 1:nSubjects
    detect = groupDataAll(iS).detectHMFC;
    discrim = groupDataAll(iS).discrimCI;
    acc = groupDataAll(iS).acc;
    rt = groupDataAll(iS).rt;
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
            groupData.rt(iV,iT,iS) = nanmean(rt{iV,iT});
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

%% P/A
for iS = 1:nSubjects
    t = behav(iS).responseTarget;
    v = behav(iS).cueValidity;
    acc = behav(iS).acc;
    
    tPA = double(behav(iS).targetType > 0);
    ntPA = double(behav(iS).nontargetType > 0);
    tPA(isnan(behav(iS).targetType)) = NaN;
    ntPA(isnan(behav(iS).nontargetType)) = NaN;
    
    PAs = [1 0]; % [pres abs]
    
    for iT = 1:numel(targets)
        wt = t==targets(iT);
        for iV = 1:numel(cueValidities)
            wv = v==cueValidities(iV);
            for iPATarget = 1:numel(PAs)
                wtpa = tPA==PAs(iPATarget);
                for iPANontarget = 1:numel(PAs)
                    wntpa = ntPA==PAs(iPANontarget);
                    
                    w = wt & wv & wtpa & wntpa;
                    
                    groupData.paAcc(iV,iT,iPATarget,iPANontarget,iS) = nanmean(acc(w));
                end
            end
        end
    end
end

%% confusion matrix
responseOptions = unique(behav(end).correctResponse(~isnan(behav(end).correctResponse)));
responseOptionsM = [0; responseOptions]; % include missed trials

confusion = nan(numel(responseOptions), numel(responseOptionsM), nSubjects);
for iSubject = 1:nSubjects
    for iP = 1:numel(responseOptions) % presented
        w = behav(iSubject).correctResponse==responseOptions(iP);
        responses = behav(iSubject).response(w);
        responses(isnan(responses)) = 0;
        for iR = 1:numel(responseOptionsM) % responded
            confusion(iP,iR,iSubject) = nnz(responses==responseOptionsM(iR))./numel(responses);
        end
    end
end

%% group summary
measures = fields(groupData);
nM = numel(measures);
for iM = 1:nM
    m = measures{iM};
    sdim = numel(size(groupData.(m))); % subject dim is always the last dim
    groupMean.(m) = mean(groupData.(m),sdim);
    groupSte.(m) = std(groupData.(m),0,sdim)./sqrt(nSubjects);
end

%% plot group data
measures = {'overallAcc','rt'};
nM = numel(measures);
figure('color','w')
for iM = 1:nM
    subplot(1,nM,iM)
    m = measures{iM};
    p1 = errorbar(groupMean.(m)', groupSte.(m)','.','MarkerSize',20);
    set(p1(2),'color','r')
    xlim([.5 2.5])
    switch m
        case 'dprimeDetect'
            ylim([0 3])
        case 'critDetect'
            ylim([-1 1])
            hold on
            plot([.5 2.5],[0 0],'--k')
        case 'rt'
            ylim([0 1.6])
        otherwise
            ylim([0 1])
    end
    title(m)
    set(gca,'XTick',[1 2])
    set(gca,'XTickLabel',{'T1','T2'})
    box off
end
legend('valid','invalid')

figure
imagesc(mean(confusion,3),[0 1])
xlabel('responded')
ylabel('presented')
set(gca,'XTick',1:numel(responseOptionsM))
set(gca,'XTickLabel',responseOptionsM)
set(gca,'YTick',1:numel(responseOptions))
colorbar

%% plot PA data
names = {'P','A'};
m = 'paAcc';
figure
for iPATarget = 1:2
    for iPANontarget = 1:2
        subplot(2,2,(iPANontarget-1)*2+iPATarget)
        hold on
        p1 = errorbar(groupMean.(m)(:,:,iPATarget,iPANontarget)', ...
            groupSte.(m)(:,:,iPATarget,iPANontarget)','.','MarkerSize',20);
        set(p1(2),'color','r')
        plot([.5 2.5],[.33 .33],'--k')
        xlim([.5 2.5])
        ylim([0 1])
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'})
        if iPATarget==1 && iPANontarget==1
            ylabel('Proportion correct')
        end
        box off
        title(sprintf('target %s, nontarget %s', names{iPATarget}, names{iPANontarget}));
%         set(gca,'LineWidth',1)
    end
end
legend('valid','invalid')
% rd_supertitle2(m);

%% plot individual data
indivM = {'overallAcc'};
nM = numel(indivM);
figure('Color','w','Position',[0 0 1200 800])
for iS = 1:nSubjects
    for iM = 1:numel(indivM)
        row = ceil(iS/2);
        if mod(iS,2)
            col = iM;
        else
%             col = iM+nM+1;
            col = iM+nM;
        end
%         subplot(ceil(nSubjects/2),nM*2+1,(row-1)*(nM*2+1)+col)
        subplot(ceil(nSubjects/2),nM*2,(row-1)*(nM*2)+col)
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
            case 'rt'
                ylim([0 1.6])
            otherwise
                ylim([0 1])
        end
        if row==1
            title(m)
        else
%             set(gca,'XTickLabel','')
%             set(gca,'YTickLabel','')
        end
%         if col==1 || col==nM+2
            ylabel(sprintf('%s',subjects{iS}))
%         end
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'})
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




