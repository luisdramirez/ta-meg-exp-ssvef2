%rd_MEG_TAPilot automatic runs
function taAutoRuns(subjectID, runs)
    % runs is a vector 
    for iRun = 1:length(runs)
        curr_run = runs(iRun);
        rd_MEG_TAPilot(subjectID, curr_run)
        sca
        % move data file movefile(source, destinataion) 
        % movefile('~/Desktop/2017...',['~/Documents/TADetectDiscrim/data/' subj]) ??
    end
    
%implement graceful exit
%auto adjust staircase after 1st run (manually attach relevent params to
%stimulus structure)

%grab all data files on desktop and move to local data for the subject

%20170318T172526_taDetectDiscrim10.mat
