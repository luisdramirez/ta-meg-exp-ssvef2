%timing irregularities
location = 'L1'; %'L1' 'laptop'
scanType = 'all frames'; % 'all frames' 'irregularity frames'



switch scanType
    case 'all frames'
        timeScanNames = {'stimulusSeq' 'targetSeq' 'targetPosSeq' 'trigSeq' 'keyCodeSeq' 'soundSeq'};
        switch location
            case 'laptop'
                timeScan = table(stimulus.seq, stimulus.target.seq, stimulus.target.posSeq, stimulus.trigSeq, stimulus.keyCodeSeq, stimulus.soundSeq);
                timeScan.Properties.VariableNames = timeScanNames;
            case 'L1'
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
            case 'L1'
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

% 
% for i = 1:length(stimulus.seq)
%     % stimulus.target.seq, 
% end
% 
% for i = 1:length(targetBlockOrder)
%     % targetBlockOrder, targetTypesBlockOrder, cueBlockOrder
% end