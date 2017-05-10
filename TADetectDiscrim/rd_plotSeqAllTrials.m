% rd_plotSeqAllTrials.m

switch p.jitter
    case 'ITI'
        blockLengths = stimulus.itiSeq + p.blockDur;
    case 'blockPrecueInterval'
        blockLengths = stimulus.jitSeq + p.blockDur;
end
blockStarts = [0 cumsum(blockLengths(1:end-1))];
blockIdx = round(blockStarts*p.refrate + 1);
nBlocks = numel(blockIdx);

precueStarts = blockStarts + .5 + stimulus.jitSeq;
precueIdx = round(precueStarts*p.refrate + 1);

% window = [0 336];
window = [-90 306];
windowSz = diff(window);
a = zeros(nBlocks,windowSz);
for i = 1:nBlocks
    win = precueIdx(i)+window(1):precueIdx(i)+window(2)-1;
    idx = 1:length(win);
    idx(win<1) = [];
    a(i,idx) = stimulus.seq(win(idx));
end

figure
subplot(3,1,1)
plot(a')
xlabel('time (frames)')
ylabel('image number in stimulus.seq')

offsets = repmat(1:nBlocks,windowSz,1);
%figure
hold on
subplot(3,1,2)
plot(a' + offsets);
xlabel('time (frames)')
ylabel('image number in stimulus.seq')