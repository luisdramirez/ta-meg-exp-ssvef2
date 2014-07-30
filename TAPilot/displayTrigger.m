function displayTrigger(trigSeq)

nTrigs = length(trigSeq);
nTrigChannels = 8;
trigMat = zeros(nTrigChannels, nTrigs);

for iTrig = 1:nTrigs
    trigVal = trigSeq(iTrig);
    
    a = dec2base(trigVal,2,nTrigChannels);
    b = nTrigChannels+1-strfind(a,'1');
    
    trigMat(b, iTrig) = 1;
end

% plot
nTimesPerPlot = 500;
nSubplots = ceil(nTrigs/nTimesPerPlot);
figure
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