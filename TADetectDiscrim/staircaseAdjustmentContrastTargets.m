function contrasts = staircaseAdjustmentContrastTargets(pedestal, contrasts, perfs)
%
% function staircaseAdjustment(contrasts, discrims)
%
% contrast is [low-contrast(decrement) high-contrast(increment)]
% discrim is [valid-discrim-targetType1(decrement)
% valid-discrim-targetType2(increment)]

%% report current values
fprintf('pedestal contrast: %.2f\n', pedestal)
fprintf('original contrasts: [%.2f %.2f]\n', contrasts(1), contrasts(2))
fprintf('performance: [%d %d]%%\n\n', round(perfs(1)*100), round(perfs(2)*100))

%% contrast decrement
% c = logspace(-2,0,30);
c = logspace(-.551,-.051,35);
[val, pIdx] = min(abs(c-pedestal)); % find c closest to pedestal
c = c(1:pIdx-1);

[val, cIdx] = min(abs(c-contrasts(1))); % find c closest to contrast
if perfs(1) > 0.90
    if cIdx < numel(c)
        cIdx = cIdx + 1; % higher is harder
    end
elseif perfs(1) < 0.65
    if cIdx > 1
        cIdx = cIdx - 1; % lower is easier
    end
%     if perfs(1) < 0.55 % extra boost for near chance performance
%         if cIdx > 1
%             cIdx = cIdx - 1; % lower is easier
%         end
%     end
else
    % do nothing
end
contrasts(1) = c(cIdx);

%% contrast increment
c = logspace(-2,0,30);
[val, pIdx] = min(abs(c-pedestal)); % find c closest to pedestal
c = c(pIdx+1:end);

[val, cIdx] = min(abs(c-contrasts(2))); % find c closest to contrast
if perfs(2) > 0.90
    if cIdx > 1
        cIdx = cIdx - 1; % lower is harder
    end
elseif perfs(2) < 0.65
    if cIdx < numel(c)
        cIdx = cIdx + 1; % higher is easier
    end
%     if perfs(2) < 0.55 % extra boost for near chance performance
%         if cIdx < numel(c)
%             cIdx = cIdx + 1; % higher is easier
%         end
%     end
else
    % do nothing
end
contrasts(2) = c(cIdx);

%% save staircase file
save('staircase.mat', 'contrasts')

%% show new values
fprintf('new contrasts: [%.2f %.2f]\n', contrasts(1), contrasts(2))


