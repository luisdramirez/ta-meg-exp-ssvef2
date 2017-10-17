% rd_analyzeTADetectDiscrimGroupBehav.m

%% setup
exptDir = pathToExpt;

% subjects = {'rd','lr','mj','af','xw'};
% startRuns = [211, 211, 221, 221, 221];

subjects = {'mj'};
startRuns = 1;

nSubjects = numel(subjects);

%% load data
for iSubject = 1:nSubjects
    sessionDir = subjects{iSubject}; 
    behavDir = sprintf('%s/analysis/%s', exptDir, sessionDir);
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
    groupMean.(m) = mean(groupData.(m),3);
    groupSte.(m) = std(groupData.(m),0,3)./sqrt(nSubjects);
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




