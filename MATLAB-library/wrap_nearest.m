% Averaging of phases with respect to nearest wrapping.
% PARAMETERS
%   phases            1-by-N or N-by-1 array of phases (in degrees).
%   degrees           boolean (Default: false)
%                     Phase is in degrees rather than radians.
%
% RETURNS
%   wrapped     1-by-N or N-by-1, the phases wrapped to minimize the
%               absolute value of the modulo-difference.
%   average     Average of the phases after wrapping to minimize the
%               absolute value of the modulo-difference.
%
% SEE ALSO
% deg2principal rad2deg unwrap meannan
%
function [phases, average] = wrap_nearest(phases, degrees)
if nargin < 2
  degrees = false;
end
if degrees
  period = 360;
else
  period = 2*pi;
end
nonnegative_range = any(phases >= period/2);

phases = mod(phases, period);
if ~nonnegative_range % initially wrap to 0 to 360 degree range instead of -180 to 180.
  which = phases > period/2;
  phases(which) = phases(which) - period;
end

wrap = [0 period -period];
for j = 1:length(phases)
  % Too simplistic, may give sub-optimal result (e.g. for [1 179 225]):
  % [~,best] = min( abs(phases(j)+wrap - phases(1)) );
  
  % Minimize variance.
  % Solves the above example and I have not seen it perform worse,
  % but not proved it to be perfect either.
  min_yet = Inf; best = 1;
  for k = 1:3
    p = phases; p(j) = p(j) + wrap(k);
    v = var(p,1);
    if v < min_yet
      best = k;
      min_yet = v;
    end
  end
  
  if j == 1 && wrap(best) ~= 0
    % Instead of wrapping first value, wrap all other values the other way
    phases(2:end) = phases(2:end) - wrap(best);
  else
    phases(j) = phases(j) + wrap(best);
  end
end

if nargout >= 2
  % Also compute average
  average = mean(phases);
end

