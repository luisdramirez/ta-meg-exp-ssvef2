function exptspath = pathToCarrascoExpts

if exist('/Volumes/purplab','dir')
    purpdir = 'purplab';
end
if exist('/Volumes/purplab-1','dir')
    purpdir = 'purplab-1';
end
if exist('/Volumes/purplab-2','dir')
    purpdir = 'purplab-2';
end

if ~exist('purpdir','var')
    error('cannot find purplab')
end

exptspath = sprintf('/Volumes/%s/EXPERIMENTS/1_Current Experiments',purpdir);