function imMasked = maskWithAnnulus(im, outerCircleSize, innerCircleSize, blurRadius, bg)

% imMasked = maskWithAnnulus(im, outerCircleSize, innerCircleSize, blurRadius, bg)
%
% Masks an image im with an annulus whose outer dimension is 
% outerCircleSize (in pixels) and inner dimension is innerCircleSize (in 
% pixels) using a background color bg. The edge of the mask is blurred by 
% blurRadius.
%
% Output is the masked image imMasked.
%
% Rachel Denison
% 5 June 2012


c1 = drawcircularblur(outerCircleSize,blurRadius);

if mod((outerCircleSize - innerCircleSize),2)
    innerCircleSize = innerCircleSize+1;
end
c2 = drawcircularblur(innerCircleSize,blurRadius);
innerCircleStartIdx = (outerCircleSize-innerCircleSize)/2;
innerCircleIdxs = innerCircleStartIdx: ...
    innerCircleStartIdx+innerCircleSize-1;
c3 = c1;
c3(innerCircleIdxs,innerCircleIdxs) = 1-c2;

c = c3;

imMasked = im.*(1-c) + bg*c;