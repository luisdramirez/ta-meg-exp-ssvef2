% stim_demo.m

displayName = 'Carrasco_L1'; % 'meg_lcd','Carrasco_L2','Carrasco_L1','Carrasco_R1'
d = loadDisplayParams(displayName);

load stimuli/taDetectDiscrim99.mat
stimulus.target.contrast = [0 .2 .8 .95];
[bg, m, t] = contrastTarget(d,stimulus.target);

% a = bg{1};
% ma = m{1,1,1}(:,:,2);
% % ami = a.*ma*.5 + a;
% ami = (a-.5).*ma*1.5 + .5;
% amd = (a-.5).*(1-ma*.8)+.5;
% figure
% subplot(1,2,1)
% imshow(ami)
% subplot(1,2,2)
% imshow(amd)

stimulus.target.pixelsPerDegree = 100;
[bg1, m1, t1] = contrastTarget(d,stimulus.target);

% Example image
im = m{1,1,4}(:,:,1); % original
im1 = m1{1,1,3}(:,:,1); % high resolution

% Display on ptb window
gray = [127 127 127];
multisample = [];
[win, rect] = Screen('OpenWindow', 0, gray, [], [], [], [], multisample);
% [win, rect] = Screen('OpenWindow', 0, gray, [0 1000 800 1600], [], [], [], multisample);

tex = Screen('MakeTexture', win, im*255);
Screen('DrawTexture', win, tex);
Screen('Flip', win);

tex = Screen('MakeTexture', win, im1*255);
Screen('DrawTexture', win, tex, [], CenterRect([0 0 size(im)], rect));
Screen('Flip', win);