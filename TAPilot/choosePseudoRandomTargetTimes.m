function targetPresentationTimes = choosePseudoRandomTargetTimes(targetPeriodDur, nTargets, targetCushion)

if nargin==0
    targetPeriodDur = 10;
    nTargets = 1;
    targetCushion = 1;
end

isi = targetPeriodDur/(nTargets+1);
possibleTargetOnsets = isi:isi:targetPeriodDur + targetCushion;
targetChoices = randperm(numel(possibleTargetOnsets));
targetOnsets = possibleTargetOnsets(targetChoices(1:nTargets));
targetOnsets = sort(targetOnsets);
minISIDiff = 0;
count = 0;
while minISIDiff < targetCushion
    jitter = normrnd(0, mean(diff(possibleTargetOnsets))*2, 1, numel(targetOnsets));
    minISIDiff = min(diff([0 sort(targetOnsets + jitter) targetPeriodDur]));
    count = count + 1;
    if count > 10000
        error('[choosePseudoRandomTargetTimes]: Max iterations (10,000) reached')
    end
end
targetPresentationTimes = sort(targetOnsets + jitter);