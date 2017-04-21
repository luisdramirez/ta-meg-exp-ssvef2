% plot performance and contrast by run

%% initial analysis
subject = 'mj';
runs = 221:230;
date = '20170421';
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
%% accuracy as a function of trial type: dec-dec, dec-inc, inc-dec, inc-inc (contrast analysis func)
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

figure
% plot valid vs invalid for dec-dec, dec-inc, inc-dec, inc-inc (8 bars
% total)
% row = validity, col = targettype , page =  pedestal, vol = target
% for each target plot
    %valid&invalid for a pedestal and both levels of that pedestal (16 bars total)
for iT = 1:2
    %for iP = 1:2
    subplot(1,2,iT)
    bar([pAcc(:,:,1,iT) pAcc(:,:,2,iT)]')
    set(gca,'XTickLabel',{'p1-dec','p1-inc','p2-dec','p2-inc'})
    if iT == 1
        ylabel('proportion correct')
    end
    ylim(ylims)
    title(sprintf('T%d',iT))
    %end
end
legend('valid','invalid')
rd_supertitle2(sprintf('%s, runs %d-%d', subject, runs(1), runs(end)))

