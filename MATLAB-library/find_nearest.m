% [indices, found] = find_nearest(matrix, value [, max_difference [, return_all_matches]])
% PARAMETERS
%   matrix  (A-by-B) numbers that constitute the haystack.
%   value  (C-by-D)  numbers to search for in the haystack.
%   max_difference   (scalar or C-by-D) only accept matches where
%                    abs(matrix(index)-value) < max_difference.
%                    (Default: 50 * max(eps(value), eps(min(min(matrix)))).)
%   return_all_matches (Default: false)
%         If false:  The returned indices has the same size as value (C-by-D)
%                    and at each position contains either an index (the closest match)
%                    or NaN (if no match within max_difference).
%         If true:   All matching incides (within max_difference) are returned
%                    and the size of indices is not fixed (F-by-1).
%                    If more than a single value was given, this means
%                    you don't always know which index belongs to which value.
% RETURNS
%   indices (C-by-D or F-by-1)
%                    Found indexes. See return_all_matches for definition.
%   found   (C-by-D or F-by-1)
%                    The value of matrix(indices), but handling
%                    the indices==NaN without causing error.
%
% SEE ALSO
%   sub2ind to convert returned index to row and column index of a matrix.
%
function [indices, found] = find_nearest(matrix, value, max_difference, return_all_matches)
if nargin < 3
  max_difference = 50 * max(eps(value), eps(min(min(abs(matrix)))));
end
if nargin < 4
  return_all_matches = false;
end

if length(max_difference) == 1
  max_difference = max_difference * ones(size(value));
elseif any(size(max_difference) ~= size(value))
  error('Different sizes for max_difference and value only allowed if max_difference is scalar');
end

% One result per value sought
if return_all_matches
  indices = [];
else
  indices = NaN(size(value));
end
for i = 1:size(value, 1)
  for j = 1:size(value, 2)
    % Search for value(i,j)
    distance = abs(matrix - value(i,j));
    which = find(distance <= max_difference(i,j) );

    if return_all_matches
      % Return all indices with distance within max_difference
      indices = [indices; which];
    elseif length(which) == 1
      % A single index found
      indices(i,j) = which;
    elseif length(which) == 0
      % None found
      indices(i,j) = NaN;
    else
      % Multiple value found, find the smallest distance (or the first among the smallest)
      [m, which] = min(distance(:));
      indices(i,j) = which;
    end

  end
end

if nargout >= 2
  % Return found data too.
  if any(any(isnan(indices)))
    % Avoid passing NaN as index
    i = indices; i(isnan(i)) = 1;
    found = matrix(i);
    found(isnan(i)) = NaN;
  else
    found = matrix(indices);
  end
end


