function coordsRot = rotateCoords(coords, rotation)
%
% function coordsRot = rotateCoords(coords, rotation)
%
% assumes the coords are centered around (0,0)
% coords have 2 rows, x and y
% rotation angle in degrees
%
% Rachel Denison
% November 2014

rotation = rotation*pi/180;

x = coords(1,:);
y = coords(2,:);

r = sqrt(x.^2+y.^2);
theta = atan2(y,x);
% [y,x] = cart2pol(theta,r)

% apply rotation
thetaRot = theta + rotation;

xRot = r.*sin(thetaRot);
yRot = r.*cos(thetaRot);
% [y,x] = pol2cart(theta,r)

coordsRot(1,:) = xRot;
coordsRot(2,:) = yRot;
