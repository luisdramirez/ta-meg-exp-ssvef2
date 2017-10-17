function exptpath = pathToExpt(directory)

% exptpath = pathToExpt(directory)

% exptpath = sprintf('%s/Luis/ta-meg-exp-ssvef2/TADetectDiscrim', pathToCarrascoExpts);
exptpath = sprintf('%s/Rachel/TA_MEG/Behav_Pilot&Training/TAContrast', pathToCarrascoExpts);

if nargin==1
    exptpath = sprintf('%s/%s', exptpath, directory);
end