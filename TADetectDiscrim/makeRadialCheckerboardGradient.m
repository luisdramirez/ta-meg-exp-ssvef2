function img = makeRadialCheckerboardGradient(pixelsPerDegree, sz, phase, contrast, thetaCycles, E, A, b, gradientSlope, gradientAngle)
%
% img = makeRadialCheckerboard(pixelsPerDegree, sz, phase, contrast, thetaCycles, E, A, b, gradientSlope, gradientAngle)

if nargin==0
    pixelsPerDegree = 100;
    sz = 5;
    thetaCycles = 8;
    E = 0.05;
    A = 1; % 0.8;
    b = 0.2; %.04;
    contrast = 0.5;
    phase = 0;
    gradientSlope = 1;
    gradientAngle = 0;
end

% Schira 2007 J Neurophys:
% M = 19.2./(abs(grid) + .77);
% plot(grid, 1./M)

% Meshgrid
grid = -sz/2:1/pixelsPerDegree:sz/2;
[x,y] = meshgrid(grid,grid);
[theta,r] = cart2pol(x,y);

% Radial checkerboard setup
checksizeR = E .* r .^ A + b;
fr = 1 ./ checksizeR;

% Polar checkerboard setup
ftheta = ones(size(theta)) .* (thetaCycles/(2*pi));

% Make checkerboard image
img = sign(sin(2 * pi * fr .* r) .* sin(2 * pi * ftheta .* theta + phase));

% Make contrast gradient
% component slopes should sum to gradientSlope
gradientAngle = mod(gradientAngle + 360,180); % deal with negative angles
yxRatio = tan(gradientAngle*pi/180);
mx = gradientSlope/(1+yxRatio);
my = gradientSlope-mx;
g = mx.*x + my.*y;

% scale images according to contrast and place within [0 1] range
img = (img * contrast) / 2 + 0.5;

