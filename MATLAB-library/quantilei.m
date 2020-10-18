% QUANTILEI Quantiles of a sample, with interpolation.
% Q = QUANTILEI(X, FREQUENCY, P) returns quantiles
% for the empirical probability distribution FREQUENCY(X).
%
% Cubic interpolation is used to give good/smooth results
% even when the distribution is sparse or noisy. The algorithm
% is the "shape-preserving piecewise cubic interpolation" of interp1.
%
% The range of Q is slightly larger than the range of the input X.
% P=0 gives X(1
%
% In intervals where the FREQUENCY is zero the cumulative distribution
% function CDF is increased by up to 5% of the minimal nonzero value
% in FREQUENCY, to handle interp1's requirement of unique CDF-values. This only
% affects the interpolation between known points, and should be negligable
% compared to the arbitraryness that sparse data (low signal-to-noise)
% itself suggests.
%
% PARAMETERS:
%   x          Coordinate to return quantile values in.        (1-by-N)
%   frequency  The frequency or probability density function,  (1-by-N)
%              each frequency-value corresponds to an x-value.
%   p          Cumulative probability values.                  (1-by-M)
%              Default: [0.03 0.1 0.25 0.5 0.75 0.9 0.97]
%   method     How the interpolation should be done. Default: 'linear'.
%              'linear' - welldefined meaning mathematically,
%              should be OK even when course sampling/noisy data.
%
%              'cubic' and 'pchip' (same) give a smooth curve but the
%              meaning (e.g. if sampled tightly and transformed to a spectrum)
%              may be a bit more meaningless & sensitive near regions of frequency==0.
%              Seems OK if the sampling is tight, so frequency does not
%              vary very much (at most 10%) between adjecent x-values.
%
%              'nearest' is not recommended, because it usually returns
%              x-coordinates halfway between 
%                 
% RETURN:
%   q          The (interpolated) x-coordinates corresponding  (1-by-M)
%              to the given cumulative probability values.
%              q(i) contains the quantile at level p(i).
%
% See also: quantile (does not handle frequency variable).
% If the appropriate form of data is given and the sampling is tight then
% quantile(xsamples,pi) =approx= round(quantilei(x,frequency,p,'linear'))
% See also: fwhmi
function q = quantilei(x, frequency, p, method)
if nargin < 3
  p = [0.03 0.1 0.25 0.5 0.75 0.9 0.97];
end
if nargin < 4
  method = 'linear';
end
  
if any(p < 0) || any(p > 1)
  error('p outside the required range from 0 to 1');
end
if any(frequency < 0)
  error('Negative probability')
end
if ~any(frequency > 0)
  disp('Warning: No positive probability')
  q = NaN(size(p));
  return;
end
if size(x,1) > 1
  error('x and frequency must be row vectors, not column vectors or matrices');
end
if any(size(x) ~= size(frequency))
  error('x and frequency must have the same size');
end
if length(x) == 1
  % When there is only one value, use it as quantlie regardless of which
  % quantile was sought
  q = x;
  return;
end

if ~issorted(x)
  [x, permutation] = sort(x);
  frequency = frequency(permutation);
  clear permutation;
end

% Prepare for calculating the cumulative probability distribution
c = cumsum(frequency);
cstep = min(frequency(frequency > 0)) / c(end);
c = c / c(end); % = c / sum(frequency)
c(c > 1) = 1; c(end) = 1; % avoid rounding errors concerning the range of c

% frequency(x) is the discrete probability function.
% (Assuming method=='linear' in the description.)
% Now, expand each sample (x,frequency) into a rectangular bar
% centered on x and reaching halfway to the previous and next x-values,
% e.g. Fbar(x) = frequency(round(x/h)*h) if x is evenly spaced by h.
% Then define the continuous cumulative distribution function as
% C(x) =def= integrate{-Inf < v <=x} Fbar(v) dv
% or equivalently, using linear interpolation,
% C(x) =linear interpolation of= {(xs,cs})
%     where
%       cs = [  0  cumsum(frequency)]
%       xs = [x_0  0.5*x_i+0.5*x_{i+1}  xs_end]
%       1 <= i < N, i integer, N = length(x)
%       x_0    = x_1 - 0.5*(x_2-x_1)    = (3*x_1 - x_2)/2
%       xs_end = x_N + 0.5*(x_N-x_{N-1}) = (3*x_N - x_{N-1})/2
%
% In the program x and c are reused also for the variables xs and cs:
x = [(3*x(1)-x(2))/2 (x(1:end-1)+x(2:end))/2  (3*x(end)-x(end-1))/2];
c = [0               c];

% Quantiles are obtained trough the inverse of the C(x) function
% mentioned above, C(x) is the interpolation of cordinates (x, c).
% However, when there are regions with frequency(x)=0 there are
% flat regions in C(x) (duplicate values in c) which means it is not
% invertible, and interp1 does not work when used for the "inverted"
% coordinates (c, x).
% The ivertability is fixed by making small slopes in the originally
% flat regions. The slope is such that the maximum change in introduced
% in c is 1% of the previously minimum step size in c (cstep).
% Additionally, duplicates at the start and the end of c are trimmed away.

% Trim away zeros if more than one, only keep one leading zero in c
first = find(c > 0, 1) - 1; %(first >= 1 because we just set c(1) = 0)
x = x(first:end);
c = c(first:end);
c0= c;%%DEBUG
% Make changes in c where necessary
d0 = find(diff(c)==0);
tic
while length(d0) > 0
  duplicate_index = d0(1);
  next_index = find(c > c(duplicate_index), 1);
  if ~isempty(next_index)
    fix_indices = (duplicate_index+1):(next_index-1);
    steps = 1:((next_index-1)-duplicate_index); % one value per index to fix
    % make the last fix be 5% of the minimal distance between c-values
    steps = steps * max(0.05*cstep/steps(end), 6*eps(c(duplicate_index+1)));;
    % IMPROVEMENT: even better would be to make steps like
    %   a*[0.1, 0.2, 0.3, ..., D-..., D-0.3, D-0.2, D-0.1]
    % where a*D = c(next_index)-c(duplicate_index) and a is
    % such that a*[0.1 0.2 0.3 ...] reproduces the current steps variable.
    % That would make the "step" in q-value be centered in the range of zero frequency.

    % make small increments in c to avoid duplicate values
    c(fix_indices) = c(fix_indices) + steps;

    df = diff(c([fix_indices next_index]));
    if any(df <= 0)
      d0 = [duplicate_index+find(df<=0) d0((next_index-duplicate_index):end)];
    else
      % It is probably enough to use the else-statement, but skipping the
      % if-statement raises the risk of infinite looping (and slowdown is negligable)
      d0 = d0((next_index-duplicate_index):end); %remove the duplicates that were fixed in this step
    end
    clear df;
  else
    % There is no next higher value, it means we have the highest x-value
    % of nonzery frequency and can trim off the tail where there are no counts.
    c = c(1:duplicate_index);
    x = x(1:duplicate_index);
    break; % done (could set d0=[];)
  end
end

% Get the quantiles
q = interp1(c, x, p, method);
% Set the end points this way too, to make sure rounding errors did not cause NaN (method='linear' does no extrapolation)
q(p <= 0) = x(1);
q(p >= 1)  = x(end);


% % A graphic comparison of interpolation methods
% pi=0:0.00005:1;
% plot(c,x,'o', pi,interp1(c,x,pi,'pchip'),'-g', p,interp1(c, x, p, 'linear'),'--r') %DEBUG
% % , pi,quantile(xxdata,pi),'.-m' 
% xlabel('Cumulative probability'); ylabel('x'); pause 

% % A test of inverting back the q-values to produce an approximative
% % probability density function (frequency), which should be like
% % the Fbar(x) described above.
% f=[18 39 83 81 6 40 53 42 66 63 0 0 0 0 0 29 43 2 99 27 21 36 40 48 56 64 68 38];
% x=1:28; pi=0:0.00008:1; xi=quantilei(x,f,pi,'linear'); %or 'cubic'
% d=diff(pi)./diff(xi); d(isnan(d) | isinf(d))=-0.01;a=sum(f);
% plot(x,f,'o-',xi,a*[0 d],'.-');
% xlabel('x'); ylabel('Frequency'); legend('Discrete f_x', 'Approximative Fbar(x)');
