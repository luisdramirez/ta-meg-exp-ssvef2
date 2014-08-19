function rd_MEG_TAPilot(n, stimfile)
% rd_MEG_TAPilot(n, stimfile)
%
% MEG SSVEP left/right flicker rate x attention expt
% with exo targets
% ------
%   Run time per experiment = 63 seconds
%       9 blocks with 7 s/block
%
% INPUTS
%   n is the runnumber [1 15]
%   stimfile is the prefix for the stimulus file containing images, stored
%       in vistadisp/Applications2/Retinotopy/standard/storedImageMatrices
%
% Modified from runme_MEG_OnOffLeftRight_ET_M2008.m
% RD, July 2014

%% Add paths
addpath(genpath('/Users/megadmin/Desktop/Experiments/Rachel/vistadisp'));

%% Settings
displayName = 'meg_lcd';
frameRate = 60;
useKbQueue = 1;
use_eyetracker = true;
eyeFile = sprintf('TA%02d%s', n, datestr(now, 'mmdd')); % 8 characters max

%% Configurations
% initialize stim tracker for MEG
PTBInitStimTracker;
global PTBTriggerLength 
PTBTriggerLength = 0.001;

% debug mode?
% PsychDebugWindowConfiguration
skipSyncTests = 0;
Screen('Preference', 'SkipSyncTests', skipSyncTests);

%% Initialize Eyetracker and do Calibration
d = loadDisplayParams('displayName',displayName,'frameRate',frameRate);
hz  = FrameRate(d.screenNumber)
if round(hz)~=frameRate
    error('Frame rate not set correctly')
end
% tr  = 1/hz*frameRate;

% Do we want to use the eyetracker?
% if n == 1; % Only for the first run
    % You have to open a screen first (to get a window pointer), to start
    % the PTBInitEyeTracker;
    d = openScreen(d);
    global PTBTheWindowPtr
    PTBTheWindowPtr = d.windowPtr;
    
    if use_eyetracker
        %Open the screen
        PTBInitEyeTracker;
        % paragraph = {'Eyetracker initialized.','Get ready to calibrate.'};
        % PTBDisplayParagraph(paragraph, {'center',30}, {'a'});
        PTBCalibrateEyeTracker;

        % actually starts the recording
        % name correponding to MEG file (can only be 8 characters!!, no extension)
        PTBStartEyeTrackerRecording(eyeFile);
    end
% end

Screen('CloseAll');

%% Default parameters
params = retCreateDefaultGUIParams;

%% Hemifield and ONOFF mixture
params.modality         = 'MEG'; 
params.calibration      = displayName;
params.experiment       = 'Experiment From File';
params.skipSyncTests    = skipSyncTests;

%% ********************
%  ***** GO ***********
%  *********************
params.loadMatrix = sprintf('%s%d.mat', stimfile, n);

% load the rest of the params, but don't start yet (rd version)
params = ret_rd(params); 

% adjust display params
% rd version has white stick on bottom to work with rotated cross
params.display = attInitFixParams_rd(params.display);
params = rotateFixCoords(params, pi/4); % rotate fix 45 deg

% set button box device
if useKbQueue
    params.display.devices.useKbQueue = 1;
    params.display.devices.keyInputInternal = [];
    params.display.devices.keyInputExternal = [];
    devices = PsychHID('devices');
    for iD=1:numel(devices)
        if strcmp(devices(iD).usageName,'Keyboard') && ...
                strcmp(devices(iD).product,'904')
            params.display.devices.keyInputExternal = iD;
        end
    end
    if isempty(params.display.devices.keyInputExternal)
        error('Did not find button box')
    end
else
    params.display.devices.useKbQueue = 0;
end

% also make sure params.devices is set the same
params.devices = params.display.devices;

% go
doRetinotopyScan(params);

%% Check timing results
f = dir('~/Desktop/2014*.mat');
fileName = fullfile('~', 'Desktop', f(end).name);
load(fileName);
figure(101); clf

% desired inter-stimulus duration
plot(diff(stimulus.seqtiming));

% measured inter-stimulus duration
hold on
plot(diff(response.flip), 'r-'); 
ylim(median(diff(response.flip)) + [-.01 .01])

% frames between stimuli
frames = round(diff(response.flip) / (1/frameRate)); 

% how many interstimulus frames differed from the median?
disp(sum(frames ~= median(frames)))

%% Check responses
figure
subplot(2,1,1)
plot(stimulus.seqtiming, stimulus.trigSeq)
title('triggers')
subplot(2,1,2)
plot(stimulus.seqtiming, response.keyCode)
title('key presses')
xlabel('seconds')

%% Rename data file with run name and number
newFileName = sprintf('%s_%s%d.mat', fileName(1:end-4), stimfile, n);
movefile(fileName, newFileName)

%% Stop Eyetracker when done with experiment
% if n == 15;
    if use_eyetracker
        PTBStopEyeTrackerRecording; % <----------- (can take a while)
    end
% end
