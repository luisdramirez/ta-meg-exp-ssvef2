% plot performance and contrast by run

%% initial analysis
subject = 'lr';
runs = 101:110;
[acc, ~, responseData_all, responseData_labels] = TADetectDiscrim_analysis(subject, runs, [], 1);
% rd_supertitle2(sprintf('%s, runs %d-%d', subject, runs(1), runs(end)))

%% mean contrast
10.^mean(log10(acc.targetContrast))

10.^mean(log10(acc.targetContrast(end/2+1:end,:)))

%% overall accuracy
% runIdx = find(strcmp(responseData_labels, 'run number'));
% correctIdx = find(strcmp(responseData_labels, 'correct'));
% nRuns = numel(unique(responseData_all(:,runIdx)));
% correct = responseData_all(:,correctIdx);
% correct(correct==-1) = 0;
% 
% for iRun = 1:nRuns
%     w = responseData_all(:,runIdx)==iRun;
%     correctByRun(iRun) = nanmean(correct(w));
% end
% correctOverall = mean(correctByRun);

% quick and dirty
accvals = acc.Discrim1_all;
correctByRun = [accvals'*[.75 .25 0 0]' accvals'*[0 0 .75 .25]'];
correctOverall = mean(correctByRun)

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
% rd_supertitle2(sprintf('%s, runs %d-%d', subject, runs(1), runs(end)))

figure
plot(acc.targetContrast)
xlabel('run')
ylabel('target contrast')

