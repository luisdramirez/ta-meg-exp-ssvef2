function img = makeRadialCheckerboardGradient(pixelsPerDegree, sz, phase, contrast, thetaCycles, E, A, b, gradientAngle, bgContrast)
%
% img = makeRadialCheckerboard(pixelsPerDegree, sz, phase, contrast, thetaCycles, E, A, b, gradientAngle, bgContrast)
%
% Schira 2007 J Neurophys:
% M = 19.2./(abs(grid) + .77);
% plot(grid, 1./M)

if nargin==0
    pixelsPerDegree = 100;
    sz = 5;
    thetaCycles = 8;
    E = 0.05;
    A = 1; % 0.8;
    b = 0.2; %.04;
    contrast = 1;
    phase = 0;
    gradientAngle = 45;
    bgContrast = 0.5;
end

unit = 'log';
gradientSlope = 1; % with respect to degrees of visual angle, later overridden by contrast

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
gradientAngle = mod(gradientAngle + 360,360); % deal with negative angles
yxRatio = tan(gradientAngle*pi/180);
mx = gradientSlope/(1+yxRatio);
my = gradientSlope-mx;
g = mx.*x + my.*y;
if gradientAngle > 135 && gradientAngle <= 315 % weird but necessary
    g = -g;
end

% Scale gradient
glims = [min(g(:)) max(g(:))];
switch unit
    case 'log'
        lbg = abs(log(bgContrast));
        g = (g/diff(glims)-0.5)*lbg*2;
        g = (g + lbg)*contrast - lbg;
        g = exp(g); % convert from log contrast units (log(c)=x, c=exp(x))        
    case 'linear'
        g = g*contrast/diff(glims);
end

% scale images according to contrast and place within [0 1] range
% img = (img .* (g + 0.5)) / 2 + 0.5;
img = (img .* g) / 2 + 0.5;

