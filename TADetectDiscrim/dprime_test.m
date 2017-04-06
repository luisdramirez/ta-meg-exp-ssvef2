% dprime_test.m

nsignaltrials = 18;
nnoisetrials = 12;

faweight = .3;

h = 0:nsignaltrials;
fa = 0:nnoisetrials;

hh = repmat(h',1,length(fa));
ff = repmat(fa,length(h),1);
ns = ones(length(h),1)*nsignaltrials;
nn = ones(length(h),1)*nnoisetrials;

for i = 1:length(fa)
    [dd(:,i), cc(:,i)] = rd_dprimeWeighted(hh(:,i), ff(:,i), ns, nn, faweight);
end

stepsizeH = diff(dd(:,1));
stepsizeFA = diff(dd(1,:));

% figure
% subplot(1,2,1)
% plot(dd)
% xlabel('n hit')
% ylabel('d-prime')
% legend(num2str(ff(1,:)'))
% subplot(1,2,2)
% plot(cc)
% xlabel('n hit')
% ylabel('c')

figure
subplot(1,2,1)
imagesc(dd)
title('d-prime')
subplot(1,2,2)
imagesc(cc)
title('c')

figure
hold on
plot(stepsizeH)
plot(-stepsizeFA')
ylabel('d-prime step size')
title(sprintf('faweight = %.1f',faweight))

a = dd>.49 & dd<.83;
figure
imagesc(a)