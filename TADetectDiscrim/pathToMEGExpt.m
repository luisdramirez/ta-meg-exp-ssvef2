function exptpath = pathToMEGExpt(directory)

% exptpath = pathToExpt(directory)

exptpath = sprintf('%s/TA_MEG/MEG/TADetectDiscrim/Behav_Data', pathToCarrascoExpts);

if nargin==1
    exptpath = sprintf('%s/%s', exptpath, directory);
end