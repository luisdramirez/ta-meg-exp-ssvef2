function f = drawcircularblur(imagesize,amt,amtin)

% function f = drawcircularblur(imagesize,amt,amtin)
%
% <imagesize> is number of pixels along one side
% <amt> is fraction of radius for blurred border on the outer edge.
%   can also be [X Y] where X is the fraction of radius offset from
%   the outer edge at which blurred border starts, and Y is the
%   fraction of radius offset from the outer edge at which blurred
%   border ends.  note that in the case that you specify a single
%   number A, this is equivalent to [A 0].
% <amtin> (optional) is [X Y] where X is the fraction of radius
%   for offset from the center and Y is the fraction of radius
%   for blurred border on the inner edge.  if [] or not supplied,
%   don't do this edge.
%
% return a matrix where values are in [0,1].

% deal with input
if length(amt)==1
  amt = [amt 0];
end
if ~exist('amtin','var') || isempty(amtin)
  amtin = [];
end

% set up grid
temp = (1:imagesize) - (1+imagesize)/2;
[xx,yy] = meshgrid(temp,temp);

% handle outer border
  % determine boundaries
bottom = imagesize/2 - (amt(1) * imagesize/2);
top = imagesize/2 - (amt(2) * imagesize/2);
  % do it
f = matrixnormalize(sqrt(xx.^2 + yy.^2),0,1,bottom,top,1);

% handle inner border
if ~isempty(amtin)
    % determine boundaries
  bottom = amtin(1) * imagesize/2;
  top = bottom + (amtin(2) * imagesize/2);
    % do it
  f = f + (1 - matrixnormalize(sqrt(xx.^2 + yy.^2),0,1,bottom,top,1));
end
