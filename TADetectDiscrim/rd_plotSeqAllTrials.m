% rd_plotSeqAllTrials.m
run = '9';
location = 'L1';
switch location
    case 'L1'
        load(['/Users/purplab/Deskptop/Luis/vistadisp/Applications2/Retinotopy/standard/storedImagesMatrices/taDetectDiscrim' run '.mat'])
end

blockLengths = stimulus.itiSeq + p.blockDur;
blockStarts = [0 cumsum(blockLengths(1:end-1))];
blockIdx = blockStarts*p.refrate + 1;
nBlocks = numel(blockIdx);

window = 250;
a = [];
for i = 1:nBlocks
    a(i,:) = stimulus.seq(blockIdx(i):blockIdx(i)+window-1);
end

%figure
subplot(3,1,1)
plot(a')
xlabel('time (frames)')
ylabel('image number in stimulus.seq')

offsets = repmat(1:nBlocks,window,1);
%figure
hold on
subplot(3,1,2)
plot(a' + offsets);
xlabel('time (frames)')
ylabel('image number in stimulus.seq')