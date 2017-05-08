% plot performance and contrast by run

%% initial analysis
subject = 'xw';
runs = 221:230;
date = '';
plotLevel = 1;
saveFile = 0;
saveFigs = 0;
[acc, ~, responseData_all, responseData_labels] = TADetectDiscrim_analysis(subject, runs, date, plotLevel, saveFile, saveFigs);
rd_supertitle2(sprintf('%s, runs %d-%d', subject, runs(1), runs(end)))

%% update behavior
b.responseData_all = responseData_all;
b.responseData_labels = responseData_labels;
b = behavior(b);

%% mean contrast
% 10.^mean(log10(acc.targetContrast))

% 10.^mean(log10(acc.targetContrast(end/2+1:end,:)))

%% overall accuracy
% quick and dirty
accvals = acc.Discrim1_all;
correctByRun = [accvals'*[.75 .25 0 0]' accvals'*[0 0 .75 .25]'];
correctOverall = mean(correctByRun);

%% accuracy as function of target type and response target
cvs = [1 -1];
for iT = 1:2
    for iTT = 1:2
        for iV = 1:2
            w = b.responseTarget==iT & b.targetType==iTT & b.cueValidity==cvs(iV);
            ttAll{iV,iTT,iT} = [b.responseTarget(w) b.targetType(w) b.cueValidity(w) b.acc(w)];
            ttAcc(iV,iTT,iT) = nanmean(b.acc(w));
        end
    end
end
%% accuracy as a function of pedestal & trial type
%pedestalSeq (below(1)/above(2) baseline), targetTypeBlockOrder (darker(1)/brighter(2)), pedestalBlockOrder
b.targetPedestal = b.responseData_all(:,end-1:end);
for iT = 1:2
    for iP = 1:2
        for iTT = 1:2
            for iV = 1:2
                w = b.responseTarget==iT & b.targetType==iTT & b.targetPedestal(:,iT)==iP & b.cueValidity==cvs(iV);
                pAll{iV,iTT,iP,iT} = [b.responseTarget(w) b.targetType(w) b.targetPedestal(w) b.cueValidity(w) b.acc(w)];
                pAcc(iV,iTT,iP,iT) = nanmean(b.acc(w)); % row = validity, col = pedestal, page = targettype, vol = target
            end
        end
    end
end

%% trial type: dec-dec, dec-inc, inc-dec, inc-inc (contrast analysis func)
%compare target type to other target type
for iT = 1:2
    for iP = 1:2
        for iNP = 1:2
            for iV = 1:2
                iNT = 3 - iT;
                w = b.responseTarget==iT & b.targetPedestal(:,iT)==iP & b.targetPedestal(:,iNT)==iNP & b.cueValidity==cvs(iV);
                typeAll{iV,iNP,iP,iT} = [b.responseTarget(w) b.targetPedestal(w) b.cueValidity(w) b.acc(w)];
                typeAcc(iV,iNP,iP,iT) = nanmean(b.acc(w));
            end
        end
    end
end

%% confusion matrix, target presented vs responded
responseOptions = unique(b.correctResponse(~isnan(b.correctResponse)));
responseOptionsM = [0; responseOptions]; % include missed trials

confusion = nan(numel(responseOptions), numel(responseOptionsM));
for iP = 1:numel(responseOptions) % presented
    w = b.correctResponse==responseOptions(iP);
    responses = b.response(w);
    responses(isnan(responses)) = 0;
    numel(responses)
    for iR = 1:numel(responseOptionsM) % responded
        confusion(iP,iR) = nnz(responses==responseOptionsM(iR))./numel(responses);
    end
end

%% confusion matrix, non-target presented vs responded
confusionNT = nan(numel(responseOptions), numel(responseOptionsM));
for iP = 1:numel(responseOptions) % presented
    w = b.swapResponse==responseOptions(iP);
    responses = b.response(w);
    responses(isnan(responses)) = 0;
    numel(responses)
    for iR = 1:numel(responseOptionsM) % responded
        confusionNT(iP,iR) = nnz(responses==responseOptionsM(iR))./numel(responses);
    end
end
        
%% plot
ylims = [0 1];
figure
subplot(1,2,1)
hold on
plot(accvals(1:2,:)')
plot(correctByRun(:,1),'k--')
ylim(ylims)
xlabel('run')
ylabel('discrim/detected')
title('T1')
subplot(1,2,2)
hold on
plot(accvals(3:4,:)')
plot(correctByRun(:,2),'k--')
ylim(ylims)
xlabel('run')
title('T2')
rd_supertitle2(sprintf('%s, runs %d-%d', subject, runs(1), runs(end)))

figure
plot(acc.targetContrast)
xlabel('run')
ylabel('target contrast')

figure
for iT = 1:2
    subplot(1,2,iT)
    bar(ttAcc(:,:,iT)')
    set(gca,'XTickLabel',{'dec','inc'})
    if iT==1
        ylabel('proportion correct')
    end
    ylim(ylims)
    title(sprintf('T%d',iT))
end
legend('valid','invalid')
rd_supertitle2(sprintf('%s, runs %d-%d', subject, runs(1), runs(end)))

figure('Position',[235 265 915 370])
for iT = 1:2
    subplot(1,2,iT)
    bar([pAcc(:,:,1,iT) pAcc(:,:,2,iT)]')
%     set(gca,'XTickLabel',{'p1-dec','p1-inc','p2-dec','p2-inc'})
    set(gca,'XTickLabel',responseOptions)
    if iT == 1
        xlabel('contrast level')
        ylabel('proportion correct')
    end
    ylim(ylims)
    title(sprintf('T%d',iT))
end
legend('valid','invalid')
rd_supertitle2(sprintf('%s, runs %d-%d', subject, runs(1), runs(end)))

figure
for iT = 1:2
    subplot(1,2,iT)
    bar([typeAcc(:,:,1,iT) typeAcc(:,:,2,iT)]')
    set(gca, 'XTickLabel',{'dec-dec','dec-inc','inc-inc','inc-dec'})
    if iT == 1
        ylabel('proportion correct')
    end
    ylim(ylims)
    title(sprintf('T%d',iT))
end
legend('valid','invalid')
rd_supertitle2(sprintf('%s, runs %d-%d', subject, runs(1), runs(end)))

figure('Position',[360 450 690 245])
subplot(1,2,1)
imagesc(confusion,[0 1])
xlabel('responded')
ylabel('presented')
title('target')
set(gca,'XTick',1:numel(responseOptionsM))
set(gca,'XTickLabel',responseOptionsM)
set(gca,'YTick',1:numel(responseOptions))
colorbar

subplot(1,2,2)
imagesc(confusionNT,[0 1])
xlabel('responded')
ylabel('presented')
title('non-target')
set(gca,'XTick',1:numel(responseOptionsM))
set(gca,'XTickLabel',responseOptionsM)
set(gca,'YTick',1:numel(responseOptions))
colorbar
