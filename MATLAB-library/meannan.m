% Mean value, ignoring NaN in each column.
% The mean of an empty array is NaN, unless the third argument emptyout=true.
% 
% PARAMETERS
%  x        M-by-N matrix of data
%  dim      the dimension to take mean along.
%           Default: the vector dimension if x is a vector,
%           or 1 if x is a  matrix.
%  emptyout (Default: false)
%           Determines what scalar mean value to return if the data is
%           empty array or all NaN.
%           If the dimensions are such that the output would be scalar
%           and emptyout==true then [] is returned, otherwise NaN is returned
%           (at the respective index if output is not scalar).
%
% RETURN
%  out   if dim==1: the means of each column (1-by-N)
%        if dim==2: the means of each row    (M-by-1)
function out = meannan(x, dim, emptyout)
if length(size(x)) > 2
  error('Higher-dimensional input than 2D matrix not supported yet.')
end
if nargin < 2
  if size(x, 1) == 1
    dim = 2; % row vector
  else
    dim = 1;
  end
end
if dim == 2
  x = transpose(x);
end
if nargin < 3
  emptyout = false;
end

M = size(x, 1);
N = size(x, 2);
out = NaN * ones(1, N);
for n = 1:N
  d = x(:,n);
  which = ~isnan(d);
  if ~any(which) && emptyout && N==1
    % Empty (or just NaN) 1D input, and desired to return empty array then
    out = [];
  else
    out(1, n) = mean(d(which));
  end
end

if dim == 2
  out = transpose(out);
end


