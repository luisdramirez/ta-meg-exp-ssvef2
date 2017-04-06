function img = makeSpiralCheckerboard(pixelsPerDegree, sz, phase, contrast, thetaCycles, E, A, b)
%
% img = makeSpiralCheckerboard(pixelsPerDegree, sz, fr, ftheta, phase, contrast)

if nargin==0
    pixelsPerDegree = 100;
    sz = 5;
%     fr = 1; % spatial frequency of a radial cycle
    thetaCycles = 8; % number of spirals
    E = 0.1;
    A = 1; 
    b = 0.2; 
    contrast = 0.5;
    phase = 0;
end

% Meshgrid
grid = -sz/2:1/pixelsPerDegree:sz/2;
[x,y] = meshgrid(grid,grid);
[theta,r] = cart2pol(x,y);

% Radial expansion (for no expansion, just set fr to a constant or E to zero)
checksizeR = E .* r .^ A + b;
fr = 1 ./ checksizeR;

% Spiral
ftheta = thetaCycles/(2*pi);
spiral = sin(2 * pi * fr .* r  +  2 * pi * ftheta .* theta + phase);

% Binarize and make checkerboard
img1 = sign(spiral);
img2 = fliplr(sign(spiral));
img = img1==img2;

% Hack for counterphase flicker
if phase==pi
    img = 1-img;
end

% Scale image according to contrast and place within [0 1] range
img = (img - 0.5) * contrast + 0.5;


