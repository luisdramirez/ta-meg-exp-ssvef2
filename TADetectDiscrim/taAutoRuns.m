%rd_MEG_TAPilot automatic runs
function taAutoRuns(runs)
    % runs is a vector 
    for iRun = 1:length(runs)
        curr_run = runs(iRun);
        rd_MEG_TAPilot(curr_run, 'taDetectDiscrim')
        sca
    end
    
%implement graceful exit
%auto adjust staircase after 1st run (manually attach relevent params to
%stimulus structure)

%grab all data files on desktop and move to local data for the subject

