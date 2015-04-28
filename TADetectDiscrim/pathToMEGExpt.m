function exptpath = pathToMEGExpt(directory)

% exptpath = pathToExpt(directory)

exptpath = sprintf('%s/TA_MEG/MEG/TADetectDiscrim', pathToCarrascoExpts);

if nargin==1
    exptpath = sprintf('%s/%s', exptpath, directory);
end