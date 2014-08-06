function f = displayTrigger(trigSeq, nBlocks)
%
% function displayTrigger(trigSeq, nBlocks)
%
% Given a trigger sequence, displays the triggers as they should appear in 
% the MEG channels
%
% trigSeq is the trigger sequence
% nBlocks (optional) is the number of blocks for plotting = number of subplots
%
% Rachel Denison
% 29 July 2014

%% deal with inputs
if nargin==1
    nBlocks = [];
end

%% make trigger channel matrix
nTrigs = length(trigSeq);
nTrigChannels = 8;
trigMat = zeros(nTrigChannels, nTrigs);

for iTrig = 1:nTrigs
    trigVal = trigSeq(iTrig);
    
    a = dec2base(trigVal,2,nTrigChannels);
    b = nTrigChannels+1-strfind(a,'1');
    
    trigMat(b, iTrig) = 1;
end

%% plot
if isempty(nBlocks)
    nTimesPerPlot = 500;
    nSubplots = ceil(nTrigs/nTimesPerPlot);
else
    nSubplots = nBlocks;
    nTimesPerPlot = round(nTrigs/nBlocks);
end

f = figure;
for i = 1:nSubplots
    subplot(nSubplots,1,i)
    if i < nSubplots
        imagesc(trigMat(:,nTimesPerPlot*(i-1)+1:nTimesPerPlot*i));
        set(gca,'XTickLabel',(50:50:nTimesPerPlot)+nTimesPerPlot*(i-1))
    else
        imagesc(trigMat(:,nTimesPerPlot*(i-1)+1:end));
        set(gca,'XTickLabel',(50:50:nTrigs-nTimesPerPlot*(i-1)) + +nTimesPerPlot*(i-1))
    end
    colormap gray
end
