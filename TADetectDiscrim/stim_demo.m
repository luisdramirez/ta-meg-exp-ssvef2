% stim_demo.m

displayName = 'Carrasco_L1'; % 'meg_lcd','Carrasco_L2','Carrasco_L1','Carrasco_R1'
d = loadDisplayParams(displayName);

load stimuli/taDetectDiscrim3.mat
stimulus.target.contrast = [.1 .9];
[bg, m, t] = contrastTarget(d,stimulus.target);

a = bg{1};
ma = m{1,1,1}(:,:,2);
ami = a.*ma*.5 + a;
% ami = (a-.5).*ma*1.5 + .5;
amd = (a-.5).*(1-ma*.8)+.5;
figure
subplot(1,2,1)
imshow(ami)
subplot(1,2,2)
imshow(amd)