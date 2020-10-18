% [w, ws] = fwhmi(x, y, fractions, kw)
% Full width at half max, or at any other fraction of the maximum.
% It also uses some smoothing and interpolation for noise stability.
%
% Linear interpolation of y is used to find the x-coordinates
% (and thus widths) for each fractional height. Then the average of ws
% is returned as w. With the default fractions (or any fractions
% symmetrically around 0.5), w is an aproximation of the FWHM.
%
% If the y values are very noisy and smoothing is used then it may happen
% that some (or all?) widths in ws are NaN, because the full_height takes
% the un-smoothed signal into account and the highest peaks of noise may have
% been reduced to less than the sought fraction for ws. If this happens to
% all fractions then w is also NaN.
% If any y-value is negative then the absolute value of the most negative
% value is added to all y-values so that all values are >= 0. However, if
% all values are initially positive then no subtraction is made.
%
% The smoothing tends to overestimate the width slightly (only tested
% on gaussians and lorentzians, in the order of half the sampling interval
% or a permille of the true fwhm for Lorentzians).
% On very smooth and densely sampled data, you can get a more accurate
% value by skipping averaging (only reading out at fraction 0.5)
% and skipping smoothing (kw = 0):  w = fwhmi(x, y, [0.5], 0)
%
% PARAMETERS
%   x   the x-coordinate in which the width is measured      (1-by-N)
%   y   the height, where half the maximum gives the width   (1-by-N)
%   fractions  the fractions of the full height at which     (1-by-M)
%       The default is [.46 .49 .50 .51 .54]. You can pass [] to use the
%       default if you only want to specify kw.
%   kw  the number of samples to use in the smoothing filter (scalar)
%       (the kernel length of the convolution). Default: 9.
%       kw will be rounded to an odd number. If kw <= 5 or the length of y
%       is less than 4*kw then no smoothing will be performed.
% RETURN
%   w    The mean width, of the widths at all the fractions  (scalar)
%   ws   The width at each fraction.                         (1-by-M)
%   xc   The center x-coordinate, taken as average between
%        * the x-coordinate of maximum y after smoothing
%        * the average among the center positions for the ranges that gave ws.
%   fxs  The left and right (early and late) x-coordinate    (2-by-M)
%        at each fractional height level. Can with some further processing
%        be used to read out risetime of a pulse (using e.g. fractions = [0.1 0.9]).
%
% See also quantilei
function [w, ws, xc, fxs] = fwhmi(x, y, fractions, kw)
% kw is the "kernel width" - determines how much smoothing to apply
if nargin < 3 || isempty(fractions)
  %fractions = [.48 .50 .52];
  fractions = [.46 .49 .50 .51 .54];
  %(DEBUG: recently fractions = [.1 .50 .9]; %was used unintentionally!)
end
if nargin < 4
  kw = 9;
else
  kw = round((kw-1)/2)*2+1; % round to odd number
end
if size(x,1) > 1 || size(y,1) > 1
  error('x and y should be row vectors');
end
if size(x,2) ~= size(y,2)
  error('x and y should have the same length.')
end
m = min(y);
if m < 0
  % Make sure all(y >= 0).
  y = y - m; % = y + abs(m);
end

if ~issorted(x)
  % Make sure (x,y) points are ordered by increasing x
  [x,perm] = sort(x);
  y = y(perm);
end
full_height = max(y);

if kw > 6 && size(y, 2) > 4*kw % only convolve when x-length is (much) greater than k-length
  % Smooth y to reduce influence of noise
  %original: k = [2 3 4*ones(1, kw-4) 3 2]; % length of kernel should be odd
  k = [.8 2.6 3.3*ones(1, (kw-5)/2) 5 3.3*ones(1, (kw-5)/2) 2.6 .8]; % length of kernel should be odd
  ys = conv(y, k/sum(k));
  ys = ys(ceil(length(k)/2):end-floor(length(k)/2));
  %DEBUG: hold off; plot([0 k 0]); pause; clf %%%
  
  % Correct for the loss at the ends where some of the convolution takes
  % zeroes from outside the original y.
  endfactors = k(1:floor(length(k)/2))  / sum(k);
  endfactors = fliplr(1 - cumsum(endfactors));
  indices = 1:length(endfactors);
  ys(indices) = ys(indices) ./ endfactors;
  ys(end+1-indices) = ys(end+1-indices) ./ endfactors;
  % Further smooth the ends by averagint in 20% of the original y there
  indices = [indices indices(end)+1];
  ys(indices) = 0.8*ys(indices) + 0.2*y(indices);
  ys(end+1-indices) = 0.8*ys(end+1-indices) + 0.2*y(end+1-indices);
  
  % Reduce the full_height a bit if y has decresed because of the smoothing.
  ys_max = max(ys);
  %full_height = 0.85*full_height + 0.15*ys_max;
  full_height = 0.5*full_height + 0.5*ys_max; % a higher fraction raises the risk that on very noisy data, no value in ys reaches 0.5*full_height
else % No smoothing can be performed
  ys = y;
  ys_max = full_height;
end

%DEBUG: hold on; plot(x, ys, '-b', x, y,'-r');% pause;

fractions = fractions * full_height;
fxs = NaN(2, length(fractions));
for i = 1:length(fractions)
  f = fractions(i);
  fi = find(ys >= f, 1, 'first');
  % Find the position at the first upwards slope crossing ys=f (or start of range if it is above)
  if fi > 1
    % The previous index (fi-1) is the last (before fi) where ys<f. Interpolate
    % to get a value for x where ys==f (may return x(fi) if ys(fi)==f).
    fxs(1,i) = x(fi-1) + (x(fi)-x(fi-1)) * (f-ys(fi-1))/(ys(fi)-ys(fi-1));
  else % the found index is 1, do not extrapolate
    if isempty(fi)
      fxs(1,i) = NaN;
    else
      fxs(1,i) = x(fi);
    end
  end
  % Find the position at the last downwards slope crossing ys=f (or end of range if it is above)
  fi = find(ys >= f, 1, 'last');
  if fi < size(ys, 2)
    % The next index (fi-1) is the frist (after fi) where ys<f. Interpolate
    % to get a value for x where ys==f (may return x(fi) if ys(fi)==f).
    fxs(2,i) = x(fi+1) + (x(fi)-x(fi+1)) * (f-ys(fi+1))/(ys(fi)-ys(fi+1));
  else % the found index is at the end, do not extrapolate
    if isempty(fi)
      fxs(2,i) = NaN;
    else
      fxs(2,i) = x(fi);
    end
  end
end

ws = diff(fxs, 1);
w = mean(ws(~isnan(ws))); % same as diff(mean(fxs, 2))

if nargout >= 3
  % Calculate a center position
  xc = mean([mean(fxs,1) mean(x(ys==ys_max))]);
end



