%timing irregularities
location = 'L1'; %'L1' 'laptop'
timeScanNames = {'sFlips' 'irrFlips' 'keyCodeSeq' 'trigSeq' 'soundSeq' 'targetSeq' 'targetPosSeq'};
msFlips = diff(response.flip);
irrFlips = find(msFlips > mean(msFlips)*1.25);
sFlips = irrFlips.* (16/1000);

switch location
    case 'laptop'
        timeScan = table(irrFlips', stimulus.keyCodeSeq(irrFlips), stimulus.trigSeq(irrFlips), stimulus.soundSeq(irrFlips), stimulus.target.seq(irrFlips), stimulus.target.posSeq(irrFlips));
        timeScan.Properties.VariableNames = timeScanNames;
    case 'L1'
        timeScan = [irrFlips', stimulus.keyCodeSeq(irrFlips), stimulus.trigSeq(irrFlips), stimulus.soundSeq(irrFlips), stimulus.target.seq(irrFlips), stimulus.target.posSeq(irrFlips)];
end

disp(timeScanNames);
disp(timeScan)