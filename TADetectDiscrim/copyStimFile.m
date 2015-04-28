function copyStimFile(run, location)

if nargin<2
    location = 'L2';
end

stimDir = 'stimuli';
stimFile = sprintf('taDetectDiscrim%d', run);
source = sprintf('%s/%s.mat', stimDir, stimFile);

switch location
    case {'L2','L1'}
        vistaDir = '~/Desktop/Rachel/vistadisp/Applications2/Retinotopy/standard/storedImagesMatrices';
    case 'laptop'
        vistaDir = '/Users/rachel/Software/vistadisp/Applications2/Retinotopy/standard/storedImagesMatrices';
    case 'MEG'
        vistaDir = '/Users/megadmin/Desktop/Experiments/Rachel/vistadisp/Applications2/Retinotopy/standard/storedImagesMatrices';
    otherwise
        error('location not recognized')
end
        
destination = sprintf('%s/%s.mat', vistaDir, stimFile);

copyfile(source, destination)