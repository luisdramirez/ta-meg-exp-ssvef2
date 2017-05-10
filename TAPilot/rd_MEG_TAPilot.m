function rd_MEG_TAPilot(run, stimfile)
commandwindow
% rd_MEG_TAPilot(run, stimfile)
% rd_MEG_TAPilot(3, 'taDetectDiscrim')

% MEG SSVEP left/right flicker rate x attention expt
% with exo targets
% ------
%   Run time per experiment = 63 seconds
%       9 blocks with 7 s/block
%
% INPUTS
%   run is the runnumber [1 15]
%   stimfile is the prefix for the stimulus file containing images, stored
%       in vistadisp/Applications2/Retinotopy/standard/storedImageMatrices
%
% Modified from runme_MEG_OnOffLeftRight_ET_M2008.m
% RD, July 2014

%% Settings
displayName = 'meg_lcd'; % 'meg_lcd', 'Carrasco_L2', 'Carrasco_L1'
frameRate = 60;
useKbQueue = 1;
use_eyetracker = false;

eyeFile = sprintf('T%02d%s', run, datestr(now, 'mmdd')); % 8 characters max
eyeDir = 'eyedata';
nStaircaseRuns = 1; % #runs for staircase to update
faWeight = 0.3;
runGUI = false; % turn GUI ON/OFF
multisample = 8;

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
% d.rect = [0 0 800 600];
d.multisample = multisample;
hz = FrameRate(d.screenNumber)
if round(hz)~=frameRate
    error('Frame rate not set correctly')
end

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
params.loadMatrix = sprintf('%s%d.mat', stimfile, run);

% load the rest of the params, but don't start yet (rd version)
params = ret_rd(params, runGUI); 

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
f = dir('~/Desktop/2017*.mat');
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

if isfield(response,'correct')
    idx = diff(response.correct)==1;
    roughAcc = response.correct(idx) + 1;
    fprintf('Estimated accuracy = %d%%\n\n', round(mean(roughAcc)*100))
end

%% Rename data file with run name and number
newFileName = sprintf('%s_%s%d.mat', fileName(1:end-4), stimfile, run);
movefile(fileName, newFileName)

%% Extra analyses if requested
if strcmp(stimfile, 'taDetectDiscrim')
    %% Analyze data from this run
    dataDir = '~/Desktop';
    vistaStimPath = 'vistadisp/Applications2/Retinotopy/standard/storedImagesMatrices';
    stimDir = sprintf('../../%s', vistaStimPath);
    plotLevel = 3; % 3 = fewest plots
    [acc, stim] = rd_analyzeTADetectDiscrimOneRun(dataDir, stimDir, run, plotLevel);

    %% Check for catch trials
    if isfield(stim.stimulus.target,'catchTrials')
        catchTrials = stim.stimulus.target.catchTrials;
    else
        catchTrials = false;
    end
    
    %% Choose the target property by which to group targets
    switch stim.p.responseOption
        case 'targetContrast4Levels'
            tg = 'targetPedestalAccAll';
            validIdx = [1 4];
            contrasts = response.target.contrast([1 4]);
            pedestal = response.target.contrast([2 3]);
        otherwise
            tg = 'targetTypeAccAll';
            validIdx = [1 3];
            contrasts = response.target.contrast;
            pedestal = stim(1).p.stimContrast;
    end
    
    %% Adjust difficulty via run-by-run staircase
    switch response.target.type
        case 'cb'
            % only use valid trials, since there's more data
            % shoot for 80% valid, mean across T1 and T2
            validDetect = mean(acc.Detect_means(validIdx));
            validDiscrim = mean(acc.Discrim1_means(validIdx));
            staircaseAdjustment(response.target.contrast, response.target.tilts(2), ...
                validDetect, validDiscrim);
        case 'contrast'
            % staircase based on 2 runs
            if nStaircaseRuns==2 && mod(run,2)==0
                % analyze previous run
                [acc(2), stim(2)] = rd_analyzeTADetectDiscrimOneRun(dataDir, stimDir, run-1, 0);
            end
            % update staircase
            if nStaircaseRuns==1 || (nStaircaseRuns==2 && mod(run,2)==0)
                for iTT = 1:2
                    validTrialsAcc{iTT} = [];
                    for iRun = 1:numel(acc)
                        for iVI = 1:numel(validIdx)
                            validTrialsAcc{iTT} = [validTrialsAcc{iTT}; acc(iRun).(tg){validIdx(iVI),iTT}];
                        end
                    end
                    validAcc(1,iTT) = nanmean(validTrialsAcc{iTT});
                end
                if catchTrials
                    validTrialsCatch = [];
                    for iRun = 1:numel(acc)
                        for iVI = 1:numel(validIdx)
                            validTrialsCatch = [validTrialsCatch; acc(iRun).catchTrialRespAll{validIdx(iVI)}];
                        end
                    end
                    staircaseAdjustmentContrastTargetsDprime(pedestal, ...
                        contrasts, validTrialsAcc, validTrialsCatch, faWeight);
                else
                    staircaseAdjustmentContrastTargets(pedestal, ...
                        contrasts, validAcc);
                end
            end
        otherwise
            error('targetType not recognized')
    end
end

%% Stop Eyetracker when done with experiment
if use_eyetracker
    PTBStopEyeTrackerRecording(eyeDir); % <----------- (can take a while)
end
