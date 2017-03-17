function f = matrixnormalize(m,minv,maxv,mmin,mmax,chop)

% function f = matrixnormalize(m,minv,maxv,mmin,mmax,chop)
%
% <m> is a matrix (no Infs).
% <minv> is the minimum desired value.  can be a scalar
%   or a matrix the same size as <m>.
% <maxv> is the maximum desired value.  can be a scalar
%   or a matrix the same size as <m>.
% <mmin> (optional) sets the min value of <m>.  can be
%   a scalar or a matrix the same size as <m>.
%   if not supplied or [], we find the actual min.
% <mmax> (optional) sets the max value of m.  can be
%   a scalar or a matrix the same size as <m>.
%   if not supplied or [], we find the actual min.
% <chop> (optional) is whether to chop off the ends
%   so that there are no values below <minv> nor
%   above <maxv>.  NaNs are preserved, however.
%   if not supplied or [], we default to 0.
%
% return <m> squeezed/stretched and moved so that [<mmin>,<mmax>] 
% maps to [<minv>,<maxv>].  if chop, we also threshold 
% non-NaN values below <minv> and values above <maxv>.
% there is a special case: if <mmin> is equal to <mmax> on 
% a per-element basis, then we just map all non-NaN values to 
% <minv>, which is an arbitrary choice.
%
% note that <chop> has no effect if <mmin> and <mmax> aren't specified.

% check empty case
if isempty(m)
  f = [];
  return;
end

% deal with input
if ~exist('mmin','var')
  mmin = [];
end
if ~exist('mmax','var')
  mmax = [];
end
if ~exist('chop','var')
  chop = 0;
end

% deal with defaults
if isempty(mmin)
  mmin = min(m(:));
end
if isempty(mmax)
  mmax = max(m(:));
end
if chop
  temp = isnan(m);
  m = max(min(m,mmax),mmin);
  m(temp) = NaN;  % preserve NaNs
end

% init f
f = zeros(size(m));

% what cases are weird?
weirdcase = mmin==mmax;
if issingle(weirdcase)
  weirdcase = repmat(weirdcase,size(m));
end

% deal with weird case
if any(weirdcase)
  f(weirdcase) = subscript(m - (mmin - minv),weirdcase);  % this just sets the values to minv, while preserving any NaNs.
end

% deal with normal case
  % want to do: f = (m-mmin) .* (maxv-minv)./(mmax-mmin) + minv
if any(~weirdcase)
  val = (maxv-minv)./(mmax-mmin);
  f(~weirdcase) = subscript(m.*val - (mmin.*val - minv),~weirdcase);  % like this for speed
end
