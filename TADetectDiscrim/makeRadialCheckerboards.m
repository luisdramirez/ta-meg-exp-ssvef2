function [img1 img2] = rd_makeRadialCheckerboards(p, r, theta)
% p = 
% r = 
% theta = 

% Radial checkerboard setup
checksizeR = p.E .* r .^ p.A;
fr = 1 ./ checksizeR;

% Polar checkerboard setup
thetaCycles = p.thetaCyc;
ftheta = ones(size(theta)) .* (thetaCycles/pi);

% Make two checkerboard images (here, contrast reversed)
img1 = sign(sin(2 * pi * fr .* r) .* sin(2 * pi * ftheta .* theta));
img2 = -sign(sin(2 * pi * fr .* r) .* sin(2 * pi * ftheta .* theta)); % cos to shift by half a cycle, sin to contrast-reverse

% scale images according to contrast and place within [0 1] range
img1 = (img1 * p.contrast) / 2 + 0.5;
img2 = (img2 * p.contrast) / 2 + 0.5;
