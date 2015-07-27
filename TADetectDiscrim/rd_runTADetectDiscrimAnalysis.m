% rd_runTADetectDiscrimAnalysis.m

subjects = {'co','ec','hl','rp','sf'};
% subjects = {'rp'};
nSubjects = numel(subjects);

runs = 411:420;

combineDates = 1; % 1 to combine all dates from each subject

plotLevel = 2; % 1 is all plots, 3 is fewest plots
saveFile = 1;
saveFigs = 1;

for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    
    switch subject
        case 'ec'
            dates = {'20150714', ...
                    '20150721', ...
                    '20150722'};
        case 'co'
            dates = {'20150716', ...
                    '20150723'};
        case 'hl'
            dates = {'20150701', ...
                    '20150702', ...
                    '20150709', ...
                    '20150710'};
        case 'rp'
            dates = {'20150625', ...
                    '20150629', ...
                    '20150702'};
%             dates = {'20150723'};
        case 'sf'
            dates = {'20150611', ...
                    '20150612', ...
                    '20150617', ...
                    '20150723'};
        otherwise
            error('subject not recognized')
    end
    
    % do analysis
    if combineDates
        TADetectDiscrim_analysis(subject, runs, dates, plotLevel, saveFile, saveFigs);
    else
        for iDate = 1:numel(dates)
            date = dates{iDate};
            TADetectDiscrim_analysis(subject, runs, date, plotLevel, saveFile, saveFigs);
        end
    end
    
    close all
end