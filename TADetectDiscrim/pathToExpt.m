function exptpath = pathToExpt(directory)

% exptpath = pathToExpt(directory)

exptpath = sprintf('%s/Luis/ta-meg-exp-ssvef2/TADetectDiscrim', pathToCarrascoExpts);

if nargin==1
    exptpath = sprintf('%s/%s', exptpath, directory);
end