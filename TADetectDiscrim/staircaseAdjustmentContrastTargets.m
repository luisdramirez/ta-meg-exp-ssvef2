function contrasts = staircaseAdjustmentContrastTargets(pedestal, contrasts, perfs)
%
% function staircaseAdjustment(contrasts, discrims)
%
% contrast is [low-contrast(decrement) high-contrast(increment)]
% discrim is [valid-discrim-targetType1(decrement)
% valid-discrim-targetType2(increment)]

%% pedestal can be 1 or 2 elements
if numel(pedestal)==1
    pedestal = [pedestal pedestal];
end

%% report current values
fprintf('pedestal contrast: [%.2f %.2f]\n', pedestal(1), pedestal(2))
fprintf('original contrasts: [%.2f %.2f]\n', contrasts(1), contrasts(2))
fprintf('performance: [%d %d]%%\n\n', round(perfs(1)*100), round(perfs(2)*100))

%% possible contrast values
% cvals = logspace(-.551,-.051,31);
cvals = logspace(-.861,-.021,22);

%% contrast decrement
c = cvals;
[val, pIdx] = min(abs(c-pedestal(1))); % find c closest to pedestal
c = c(1:pIdx-1);

[val, cIdx] = min(abs(c-contrasts(1))); % find c closest to contrast
if perfs(1) > 0.85 %0.90
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
c = cvals;
[val, pIdx] = min(abs(c-pedestal(2))); % find c closest to pedestal
c = c(pIdx+1:end);

[val, cIdx] = min(abs(c-contrasts(2))); % find c closest to contrast
if perfs(2) > 0.85 %0.90
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

%% show new values
fprintf('new contrasts: [%.2f %.2f]\n', contrasts(1), contrasts(2))


%% if 4 contrast levels
if pedestal(1)<pedestal(2)
    contrasts = [contrasts(1) pedestal contrasts(2)];
end

%% save staircase file
save('staircase.mat', 'contrasts')

