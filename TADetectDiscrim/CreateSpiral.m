function spiral = CreateSpiral(scr, gaborsiz, fre, phase, contrast)
%%
% scr: Screen parameters
% gaborfre: gabor width in visual angle
% gaborstd: std of the Gaussian modulator
% fre: gabor frequency: cycle per degree
% phase: phase
% contrast: contrast of the gabor
% There is no orientation parameter here. Rotate the texture when put on the scree.
%%
visiblefre=angle2pix(scr, [gaborsiz gaborsiz]); %angle2pix function?
%[x,y]=meshgrid(-1*visiblefre/2:visiblefre/2, -1*visiblefre/2:1*visiblefre/2);
[x,y]=meshgrid(1:visiblefre, 1:visiblefre);
x = x-mean(x(:));
y = y-mean(y(:));
gaborsiz = pix2angle(scr,visiblefre(1));
x = Scale(x)*gaborsiz-gaborsiz/2;
y = Scale(y)*gaborsiz-gaborsiz/2;

[th,r] = cart2pol(x,y);

%modulator    =exp(-((x/gaborstd).^2)-((y/gaborstd).^2));
spiral        = cos(r*fre*2*pi+phase).*contrast;

return

scr.dist = 60; %cm
scr.width = 44.5; %cm
scr.resolution = [1680,1050];
size=2;
fre=3;
gaborstd=.18
ori=0;
phase=0;
contrast=1;
spiral = CreateSpiral(scr, gaborsiz, fre, phase, contrast);
