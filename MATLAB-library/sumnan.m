% Compute sum, ignoring NaN-values.
% 
% PARAMETERS
%  x        M-by-N matrix of data
%  dim      the dimension to take mean along.
%           Default: the vector dimension if x is a vector,
%           or 1 if x is a  matrix.
%  NaN_if_none (Default: false)
%           If there are no non-NaN values to sum, should NaN or zero be
%           returned? The default (NaN_if_none==false) means zero.
% RETURN
%  out   if dim==1: the sum of each column (1-by-N)
%        if dim==2: the sum of each row    (M-by-1)
function out = sumnan(x, dim, NaN_if_none)
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
  NaN_if_none = false;
end

M = size(x, 1);
N = size(x, 2);
out = NaN * ones(1, N);
for n = 1:N
  d = x(:,n);
  which = ~isnan(d);
  if ~any(which)
    if NaN_if_none
      % Empty (or just NaN) 1D input, and desired to return empty array then
      out(1, n) = NaN;
    else
      out(1, n) = 0;
    end
  else
    out(1, n) = sum(d(which));
  end
end

if dim == 2
  out = transpose(out);
end


