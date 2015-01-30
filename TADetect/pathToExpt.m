function exptpath = pathToExpt(directory)

% exptpath = pathToExpt(directory)

exptpath = sprintf('%s/TA_MEG/TADetect', pathToCarrascoExpts);

if nargin==1
    exptpath = sprintf('%s/%s', exptpath, directory);
end