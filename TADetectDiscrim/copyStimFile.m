function copyStimFile(run)

stimDir = 'stimuli';
stimFile = sprintf('taDetectDiscrim%d', run);
source = sprintf('%s/%s.mat', stimDir, stimFile);

vistaDir = '/Local/Users/purplab/Desktop/Rachel/vistadisp/Applications2/Retinotopy/standard/storedImagesMatrices';
destination = sprintf('%s/%s.mat', vistaDir, stimFile);

copyfile(source, destination)