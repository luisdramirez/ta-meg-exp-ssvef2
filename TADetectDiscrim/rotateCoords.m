function coordsRot = rotateCoords(coords, rotation)
%
% function coordsRot = rotateCoords(coords, rotation)
%
% assumes the coords are centered around (0,0)
% coords have 2 rows, x and y
% rotation angle in degrees, positive is clockwise
%
% Rachel Denison
% November 2014

rotation = rotation*pi/180;

x = coords(1,:);
y = coords(2,:);

r = sqrt(x.^2+y.^2);
theta = atan2(y,x);

% apply rotation
% unit circle goes counterclockwise, but we want clockwise rotation, so
% subtract
thetaRot = theta - rotation; 

xRot = r.*cos(thetaRot);
yRot = r.*sin(thetaRot);
% [y,x] = pol2cart(theta,r)

coordsRot(1,:) = xRot;
coordsRot(2,:) = -yRot; % y is negative in PTB
