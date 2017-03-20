function contrasts = staircaseAdjustmentContrastTargets(pedestal, contrasts, discrims)
%
% function staircaseAdjustment(contrasts, discrims)
%
% contrast is [low-contrast(decrement) high-contrast(increment)]
% discrim is [valid-discrim-targetType1(decrement)
% valid-discrim-targetType2(increment)]

%% report current values
fprintf('pedestal contrast: %.2f\n', pedestal)
fprintf('original contrasts: [%.2f %.2f]\n', contrasts(1), contrasts(2))
fprintf('discrimination: [%d %d]%%\n\n', round(discrims(1)*100), round(discrims(2)*100))

%% contrast decrement
c = logspace(-2,0,30);
[val, pIdx] = min(abs(c-pedestal)); % find c closest to pedestal
c = c(1:pIdx-1);

[val, cIdx] = min(abs(c-contrasts(1))); % find c closest to contrast
if discrims(1) > 0.90
    if cIdx < numel(c)
        cIdx = cIdx + 1; % higher is harder
    end
elseif discrims(1) < 0.55
    if cIdx > 1
        cIdx = cIdx - 1; % lower is easier
    end
%     if discrims(1) < 0.55 % extra boost for near chance performance
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
if discrims(2) > 0.90
    if cIdx > 1
        cIdx = cIdx - 1; % lower is harder
    end
elseif discrims(2) < 0.55
    if cIdx < numel(c)
        cIdx = cIdx + 1; % higher is easier
    end
%     if discrims(2) < 0.55 % extra boost for near chance performance
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


