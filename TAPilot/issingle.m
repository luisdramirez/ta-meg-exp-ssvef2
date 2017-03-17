function f = issingle(x)

% function f = issingle(x)
%
% <x> is a matrix
%
% return whether <x> is not a struct,
% not a cell matrix, and has dimension [1 1].
%
% note that issingle([]) is 0.

f = ~isstruct(x) && ~iscell(x) && isequal(size(x),[1 1]);
