function f = subscript(m,range)

% function f = subscript(m,range)
%
% <m> is a matrix or a string referring to a matrix
%   variable in the base workspace.
% <range> is:
%   (1) vector index range that does not use
%       the 'end' keyword, e.g., 4, [5 6]
%   (2) the string ':'
%   (3) a cell vector, where elements can be ':',
%       e.g., {':' 4:5}
%
% the idea is to simply return something like m(range).
%
% this function is useful for cases where you want
% to get a range of elements from something that
% already has parentheses or brackets.
% it is also useful for working with string-variable
% representations of matrices.  it is also useful
% for using cell vector indices on-the-fly.
%
% see also flatten.m, subscriptc.m

if iscell(range)
  if ischar(m)
    f = evalin('base',['subscript(',m,',',cell2str(range),');']);
  else
    f = m(range{:});
  end
else
  if ischar(m)
    if isequal(range,':') && ischar(range)
      f = evalin('base',[m,'(:);']);
    else
      f = evalin('base',[m,'(',mat2str(range),');']);
    end
  else
    if isequal(range,':') && ischar(range)
      f = m(:);
    else
      f = m(range);
    end
  end
end
