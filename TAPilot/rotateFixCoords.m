function params = rotateFixCoords(params, rotation)

% rotation = pi/4;

fixCoords = params.display.fixCoords;
center = params.display.numPixels/2;

for i = 1:numel(fixCoords)
    coords = fixCoords{i};
    nCoords = size(coords,2);
    coordsRot = zeros(size(coords));
    
    coords0 = coords - repmat(center',1,nCoords);
    x = coords0(1,:);
    y = coords0(2,:);
    
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
    
    fixCoordsRot{i} = coordsRot + repmat(center',1,nCoords);
end

params.display.fixCoords = fixCoordsRot;