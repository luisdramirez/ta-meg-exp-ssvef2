% imageCheck.m

% load stimulus

blockOrder = order.blockOrder;
noise = p.noise;
stimSize = p.stimSize;
contrast = p.contrasts(1);
orientation = p.orientation;
spatialFreq = p.spatialFreq;
pixelsPerDegree = p.pixelsPerDegree;
blurRadius = p.blurRadius;
backgroundColor = p.backgroundColor;

iContrast = 1;
for iPhase = 1:2
    for iTrial = 1:length(blockOrder)
        luminanceDiff = 1;
        while luminanceDiff > noise.luminanceDiffThresh
            im1 = makeFilteredNoise2(stimSize, contrast/2, orientation, ...
                noise.orientationBandwidth, spatialFreq/noise.spatialFreqBandwidthFactor, ...
                spatialFreq*noise.spatialFreqBandwidthFactor, pixelsPerDegree, 0);
            im2 = makeFilteredNoise2(stimSize, contrast/2, orientation+90, ...
                noise.orientationBandwidth, spatialFreq/noise.spatialFreqBandwidthFactor, ...
                spatialFreq*noise.spatialFreqBandwidthFactor, pixelsPerDegree, 0);
            im = (im1 - 0.5) + (im2 - 0.5) + 0.5;
            immasked = maskWithAnnulus(im, length(im), 0, blurRadius, backgroundColor); % only for luminance testing
            luminanceDiff = abs(mean(immasked(:))-.5);
        end
        if iPhase==1
            s{iPhase, iContrast, iTrial} = im;
        elseif iPhase==2
            s{iPhase, iContrast, iTrial} = 1 - s{1, iContrast, iTrial};
        end
    end
end
ims = maskWithAnnulus(s{1,1,1}, length(im), 0, blurRadius, backgroundColor);
bgim = ims(:,:,1)*255;
targ0 = buildColorGrating(target.pixelsPerDegree, [target.imSize target.imSize], ...
    target.spatialFreq, tilt, phase, target.contrast, 0, 'bw');
if size(bgim,1)~=size(targ0,1)
    if size(targ0,1)==size(bgim,1)+2
        targ0 = targ0(2:end-1,2:end-1);
    else
        targ0 = targ0(2:end,2:end);
    end
end
targ1 = (targ0-.5) + (bgim/255-.5) + .5;
targ = maskWithAnnulus(targ1, size(bgim,1), 0, ...
    target.blurRadius, target.backgroundColor);
targs = targ; 

% if images already made, skip to here
ims = stimulus.images/255;
target = stimulus.target;
target.contrast = .4;
tilts = [-5 5];

tt = target.targetTypes;
tp = target.targetPedestals;
ntrials = length(stimulus.itiSeq);
blankidx = 1:5:ntrials;
targidx = setdiff(1:ntrials, blankidx);
for iT = 1:numel(tt)
    trialNum = ceil(iT/2);
    tilt = (tp(iT)-1)*90 + tilts(tt(iT));
    for iPhase = 1:2
        phase = target.phases(iPhase);
        targ0 = buildColorGrating(target.pixelsPerDegree, [target.imSize target.imSize], ...
            target.spatialFreq, tilt, phase, target.contrast, 0, 'bw');
        bgim = stimulus.images(:,:,targidx(trialNum)*2-1+iPhase-1);
        if size(bgim,1)~=size(targ0,1)
            if size(targ0,1)==size(bgim,1)+2
                targ0 = targ0(2:end-1,2:end-1);
            else
                targ0 = targ0(2:end,2:end);
            end
        end
        targ1 = (targ0-.5) + (bgim/255-.5) + .5;
        targ = maskWithAnnulus(targ1, size(bgim,1), 0, ...
            target.blurRadius, target.backgroundColor);
        targs(:,:,iT) = targ; %%% debug
        % target.textures(iT,iPhase) = Screen('MakeTexture', display.windowPtr, targ.*255);
    end
end



for i = 1:size(ims,3)
    im = ims(:,:,i);
    imContrast(i) = std(im(:));
    imMax(i) = max(im(:));
    imMin(i) = min(im(:));
    imFFT(:,:,i) = fftshift(fft2(im));
end

for i = 1:size(targs,3)
    t = targs(:,:,i);
    tContrast(i) = std(t(:));
    tMax(i) = max(t(:));
    tMin(i) = min(t(:));
    tFFT(:,:,i) = fftshift(fft2(t));
end

fs = p.pixelsPerDegree;
N = size(im,1);
df = fs/N;
sampleIndex = -N/2:N/2-1;
f = sampleIndex*df;

im = imFFT(:,:,1);
% im = fftshift(fft2(im0));
figure
imagesc(abs(im));
set(gca,'Clim',[0 80])
colorbar
set(gca,'XTick',11:10:size(im,1))
set(gca,'XTickLabel',round(f(11:10:size(im,1))*10)/10)
set(gca,'YTick',11:10:size(im,1))
set(gca,'YTickLabel',round(f(11:10:size(im,1))*10)/10)
xlabel('spatial frequency (cpd)')
% xtickformat('%.1f')

nbins = 20;
[Zr, R] = radialavg(abs(im),nbins,0,0);
figure
plot(R(2:end),Zr(2:end))
set(gca,'XTickLabel',round((0:f(end)/5:f(end))*10)/10)
xlabel('spatial frequency (cpd)')
ylabel('amplitude')


fH = figure;
subplot(4,1,1)
histogram(imContrast)
ylabel('rms')
title('noise images')
subplot(4,1,2)
histogram((imMax-imMin)./(imMax+imMin))
ylabel('michelson')
subplot(4,1,3)
histogram(imMax)
ylabel('max')
subplot(4,1,4)
histogram(imMin)
ylabel('min')
for iAx = 1:numel(fH.Children)
    fH.Children(iAx).XLim = [0 1];
    fH.Children(iAx).Children.BinWidth = .02;
end

fH = figure;
subplot(4,1,1)
histogram(tContrast)
ylabel('rms')
title('target + noise')
subplot(4,1,2)
histogram((tMax-tMin)./(tMax+tMin))
ylabel('michelson')
subplot(4,1,3)
histogram(tMax)
ylabel('max')
subplot(4,1,4)
histogram(tMin)
ylabel('min')
for iAx = 1:numel(fH.Children)
    fH.Children(iAx).XLim = [0 1];
    fH.Children(iAx).Children.BinWidth = .02;
end




