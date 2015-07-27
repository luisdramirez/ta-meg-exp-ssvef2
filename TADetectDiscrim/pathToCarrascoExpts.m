function exptspath = pathToCarrascoExpts

if exist('/Volumes/purplab','dir')
    purpdir = 'purplab';
elseif exist('/Volumes/purplab-1','dir')
    purpdir = 'purplab-1';
else
    error('cannot find purplab')
end

exptspath = sprintf('/Volumes/%s/EXPERIMENTS/1_Current Experiments/Rachel',purpdir);