function contrasts = staircaseAdjustmentContrastTargetsDprime(pedestal, contrasts, trialacc, trialcatch, faweight)
%
% function staircaseAdjustment(contrasts, discrims)
%
% contrast is [low-contrast(decrement) high-contrast(increment)]
% discrim is [valid-discrim-targetType1(decrement)
% valid-discrim-targetType2(increment)]

%% inputs
if nargin<5
    faweight = 1;
end

%% performance metrics
acc = [mean(trialacc{1}) mean(trialacc{2})];
tt1report = mean(trialcatch==1);
for iTT = 1:numel(trialacc)
    nh = sum(trialacc{iTT});
    nsignal = numel(trialacc{iTT});
    nfa = sum(trialcatch==iTT);
    nnoise = numel(trialcatch);
    [dprime(iTT), criterion(iTT)] = rd_dprimeWeighted(nh, nfa, nsignal, nnoise, faweight);
end
perfs = dprime;

%% report current values
fprintf('pedestal contrast: %.2f\n', pedestal)
fprintf('original contrasts: [%.2f %.2f]\n', contrasts(1), contrasts(2))
fprintf('valid trial accuracy: [%d %d]%%\n', round(acc(1)*100), round(acc(2)*100))
fprintf('target type 1 bias: %.2f\n', tt1report)
fprintf('valid trial d-prime: [%1.2f %1.2f]\n\n', dprime(1), dprime(2))

%% possible contrast values
cvals = logspace(-.551,-.051,31);

%% contrast decrement
c = cvals;
[val, pIdx] = min(abs(c-pedestal)); % find c closest to pedestal
c = c(1:pIdx-1);

[val, cIdx] = min(abs(c-contrasts(1))); % find c closest to contrast
if perfs(1) > 0.83
    if cIdx < numel(c)
        cIdx = cIdx + 1; % higher is harder
    end
elseif perfs(1) < 0.49
    if cIdx > 1
        cIdx = cIdx - 1; % lower is easier
    end
else
    % do nothing
end
contrasts(1) = c(cIdx);

%% contrast increment
c = cvals;
[val, pIdx] = min(abs(c-pedestal)); % find c closest to pedestal
c = c(pIdx+1:end);

[val, cIdx] = min(abs(c-contrasts(2))); % find c closest to contrast
if perfs(2) > 0.83
    if cIdx > 1
        cIdx = cIdx - 1; % lower is harder
    end
elseif perfs(2) < 0.49
    if cIdx < numel(c)
        cIdx = cIdx + 1; % higher is easier
    end
else
    % do nothing
end
contrasts(2) = c(cIdx);

%% save staircase file
save('staircase.mat', 'contrasts')

%% show new values
fprintf('new contrasts: [%.2f %.2f]\n', contrasts(1), contrasts(2))


