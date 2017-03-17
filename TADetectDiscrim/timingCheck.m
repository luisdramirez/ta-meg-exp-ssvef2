%timing irregularities

timeScanNames = {'irrFlips' 'keyCodeSeq' 'trigSeq' 'soundSeq' 'targetSeq' 'targetPosSeq'};
msFlips = diff(response.flip);
irr = find(msFlips > mean(msFlips)*1.25);

timeScan = table(irr', stimulus.keyCodeSeq(irr), stimulus.trigSeq(irr), stimulus.soundSeq(irr), stimulus.target.seq(irr), stimulus.target.posSeq(irr));
timeScan.Properties.VariableNames = timeScanNames;
disp(timeScan)