%rd_MEG_TAPilot automatic runs
function taAutoRuns(runs)
    % runs is a vector 
    for iRun = 1:length(runs)
        curr_run = runs(iRun);
        rd_MEG_TAPilot(curr_run, 'taDetectDiscrim')
    end
