function staircaseAdjustment(contrast, tilt, detect, discrim)
%
% function staircaseAdjustment(contrast, tilt, detect, discrim)

%% report current values
fprintf('original contrast: %.2f\n', contrast)
fprintf('original tilt: %1.1f\n\n', tilt)

fprintf('detection: %d%%\n', round(detect*100))
fprintf('discrimination: %d%%\n\n', round(discrim*100))

%% adjust contrast
updateTilt = 0; % only update tilt if contrast is in a good range

c = round(logspace(-1,0,10)*100)/100;
[val, cIdx] = min(abs(c-contrast)); % find c closest to contrast
if detect > 0.90
    if cIdx > 1
        cIdx = cIdx - 1;
    end
elseif detect < 0.75
    if cIdx < numel(c)
        cIdx = cIdx + 1;
    end
else
    updateTilt = 1;
end
contrast = c(cIdx);

%% adjust tilt
if updateTilt
    if tilt <= 2
        if discrim > .90
            tilt=tilt-.2;
        elseif discrim > .85
            tilt=tilt-.1;
        elseif discrim < .70
            tilt=tilt+.2;
        elseif discrim < .75 && discrim >= .70
            tilt=tilt+.1;
        else
            tilt = tilt;
        end
    elseif tilt > 2 && tilt <= 5
        if discrim > .90
            tilt=tilt-1;
        elseif discrim > .85
            tilt=tilt-.5;
        elseif discrim < .70
            tilt=tilt+1;
        elseif discrim < .75 && discrim >= .70
            tilt=tilt+.5;
        else
            tilt = tilt;
        end
    elseif tilt > 5
        if discrim > .90
            tilt=tilt-2;
        elseif discrim > .85
            tilt=tilt-1;
        elseif discrim < .70
            tilt=tilt+1.5;
        elseif discrim < .75 && discrim >= .70
            tilt=tilt+1;
        else
            tilt = tilt;
        end
    end
end

%% save staircase file
save('staircase.mat', 'contrast', 'tilt')

%% show new values
fprintf('new contrast: %.2f\n', contrast)
fprintf('new tilt: %1.1f\n\n', tilt)


