function exptpath = pathToMEGExpt(directory)

% exptpath = pathToExpt(directory)

exptpath = sprintf('%s/Rachel/TA_MEG/MEG/TANoise/Behav_Data', pathToCarrascoExpts);
% exptpath = sprintf('%s/Rachel/TA_MEG/MEG/TAContrast/Behav_Data', pathToCarrascoExpts);
% exptpath = sprintf('%s/TA_MEG/MEG/TADetectDiscrim/Behav_Data', pathToCarrascoExpts);

if nargin==1
    exptpath = sprintf('%s/%s', exptpath, directory);
end