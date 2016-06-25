% rd_runTADetectDiscrimAnalysis.m

% subjects = {'co','ec','hl','rp','sf','mr','jz','ps','nms','yh','jp','pv','nc','rl'};
subjects = {'rl'};
nSubjects = numel(subjects);

runs = 411:419;

combineDates = 0; % 1 to combine all dates from each subject

plotLevel = 1; % 1 is all plots, 3 is fewest plots
saveFile = 0;
saveFigs = 1;

for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    
    switch subject
        case 'sf'
            dates = {'20150611', ...
                '20150612', ...
                '20150617', ...
                '20150723'};
        case 'hl'
            dates = {'20150701', ...
                '20150702', ...
                '20150709', ...
                '20150710'};
        case 'ec'
            dates = {'20150714', ...
                '20150721', ...
                '20150722', ...
                '20150804'};
        case 'co'
            dates = {'20150716', ...
                '20150723', ...
                '20150811'};
        case 'rp'
            dates = {'20150625', ...
                '20150629', ...
                '20150702'};
            % dates = {'20150723'};
        case 'mr'
            dates = {'20150730', ...
                '20150804', ...
                '20150806'};
        case 'jz'
            dates = {'20150820', ...
                '20150825', ...
                '20150827'};
        case 'ps'
            dates = {'20150825', ...
                '20150827', ...
                '20150828'};
        case 'nms'
            dates = {'2015082*', ...
                '20150831', ...
                '20150902'};
        case 'yh'
            dates = {'20151112'};
        case 'yz'
            dates = {'20151112', ...
                '20151113'};
        case 'jp'
            dates = {'20151116'};
        case 'rl'
%             dates = {'20151201', ...
%                 '20151202', ...
%                 '20151210'};
%             dates = {'20151201'}; % 411-420
%             dates = {'20151202'}; % 411-418
            dates = {'20151210'}; % 411-419
        case 'pv'
            dates = {'20151203', ...
                '20151207', ...
                '20151209'};
        case 'nc'
            dates = {'20151207', ...
                '20151208'};
        otherwise
            error('subject not recognized')
    end
    
    % do analysis
    if combineDates
        TADetectDiscrim_analysis(subject, runs, dates, plotLevel, saveFile, saveFigs);
%         close all
    else
        for iDate = 1:numel(dates)
            date = dates{iDate};
            TADetectDiscrim_analysis(subject, runs, date, plotLevel, saveFile, saveFigs);
%             close all
        end
    end
end