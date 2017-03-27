%timing irregularities
location = 'L1'; %'L1' 'laptop'
scanType = 'all frames'; % 'all frames' 'irregularity frames'


%% Plot timing
switch scanType
    case 'all frames'
        timeScanNames = {'stimulusSeq' 'targetSeq' 'targetPosSeq' 'trigSeq' 'keyCodeSeq' 'soundSeq'};
        switch location
            case 'laptop'
                timeScan = table(stimulus.seq, stimulus.target.seq, stimulus.target.posSeq, stimulus.trigSeq, stimulus.keyCodeSeq, stimulus.soundSeq);
                timeScan.Properties.VariableNames = timeScanNames;
            otherwise
                timeScan = [stimulus.seq stimulus.target.seq stimulus.target.posSeq stimulus.trigSeq stimulus.keyCodeSeq stimulus.soundSeq];
        end
    case 'irregularity frames'
        timeScanNames = {'sFlips' 'irrFlips' 'keyCodeSeq' 'trigSeq' 'soundSeq' 'targetSeq' 'targetPosSeq'};
        msFlips = diff(response.flip);
        irrFlips = find(msFlips > mean(msFlips)*1.25); %location of flips that are more than 16ms
        sFlips = irrFlips.* (16/1000); %where in seconds these irregular flips happen
        switch location
            case 'laptop'
                timeScan = table(sFlips', irrFlips', stimulus.keyCodeSeq(irrFlips), stimulus.trigSeq(irrFlips), stimulus.soundSeq(irrFlips), stimulus.target.seq(irrFlips), stimulus.target.posSeq(irrFlips));
                timeScan.Properties.VariableNames = timeScanNames;
            otherwise
                timeScan = [sFlips', irrFlips', stimulus.keyCodeSeq(irrFlips), stimulus.trigSeq(irrFlips), stimulus.soundSeq(irrFlips), stimulus.target.seq(irrFlips), stimulus.target.posSeq(irrFlips)];
        end
end

% disp(timeScanNames)
% disp(timeScan)

for i=1:size(timeScan,2)  
    subplot(size(timeScan,2),1,i)
    plot(timeScan(:,i))
    title(timeScanNames{i})
    hold on
end
xlabel('frames')

%% Verify each frame 

seconds = stimulus.seqtiming;

%plots for trigger
trigFrameIndx = find(stimulus.trigSeq > 0); %when trigger is sent
trigFrameSecs = trigFrameIndx*(1/60);

soundFrameIndx = find(stimulus.soundSeq > 0); %when cues play
soundFrameSecs = soundFrameIndx*(1/60);

keyFrameIndx = find(stimulus.keyCodeSeq > 0); %when key code
keyFrameSecs = keyFrameIndx*(1/60);
responseDuration = diff(keyFrameIndx);
responseDurationSec = responseDuration*1/60;


targFrameIndx = find(stimulus.seq == 3 | stimulus.seq == 4); % target present (phase of target)
condFrameIndx = find(stimulus.target.seq == 1 | stimulus.target.seq == 2); %target present (condition of target)
posFrameIndx = find(stimulus.target.posSeq > 0);

% timing for target displays 
targFrames = [diff(targFrameIndx) diff(condFrameIndx) diff(posFrameIndx)];
targFramesSec = targFrames.*(1/60);
targFramesInc = 0;

% for i=1:length(soundFrameIndx)
%     soundTrigFrameIndx(i) = find(trigFrameIndx == soundFrameIndx(i));
% end
% 
% for i=1:length(trigFrameIndx)
%     targTrigFrameIndx(i) = find(trigFrameIndx(i) == targFrameIndx);
% end
flipFrames = diff(response.flip);
badFramesIndx = find(flipFrames > median(flipFrames)*2);
badFrames = flipFrames(badFramesIndx);

% plot(stimulus.trigSeq,'b')

plot(stimulus.soundSeq*100,'m')
hold on
plot(stimulus.target.seq*100,'g')
plot(stimulus.keyCodeSeq,'k')
plot(diff(response.flip)*1000,'r')
legend('soundSeq','targSeq','keyCodeSeq','flips')


%% Verify each trial

%targetBlockOrder, cueBlockOrder, targetTypeBlockOrder
    %targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres','abs-abs'};
    %cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2'}; % 2-1 = cueT2,postcueT1
    %targetTypeBlockOrder = 2 x 41, T1 T2 decrement (1) / increment(2)
% 
% for trial = 1:length(order.trialBlockOrder)
%     
% end