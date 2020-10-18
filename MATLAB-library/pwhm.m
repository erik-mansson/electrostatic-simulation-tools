%
%
% PARAMETERS
%  time               N-by-1 the times at which raw_signal is given [scan units].
%  raw_signal         N-by-1 the data to fit an oscillation to.
%  period             The period of the oscillation [scan units].
%  seleted_phase      The phase [degrees] where amplitude is maximum, or is preferred by other kind of optimization.
% RETURNS:
%  fullwidth_phase    Phase range within which amplitude is >= half of the amplitude for selected_phase
%  amp                (1-by-P)
%  phase              (1-by-P)
%  const              (1-by-P)
% EXAMPLE:
%  t=0:0.1:1; s = 5 + 2*cos(2*pi*t + 2*pi/180); pwhm(t, s, 1, 0)
%
function [fullwidth_phase, peak_phase, all] = pwhm(time, raw_signal, period, selected_phase)

if size(raw_signal,1) == 1 && size(raw_signal,2) > 1
  % raw_signal should be a column vector
  raw_signal = transpose(raw_signal);
end
if size(time,1) == 1 && size(time,2) > 1
  % time should be a column vector
  time = transpose(time);
end

fullwidth_phase = NaN;
N = length(raw_signal);

% Slow algorithm: evaluate all phases, so also minimum amplitude approximated. (a faster would be to start at selected_phase and step/binsearch outwards)
search_step = 1; %[degrees]
search_step = 15; % DEBUG
phase = ((-180+search_step):search_step:180) + selected_phase;
%phase = -20:10:20; % DEBUG
%disp(mat2str(phase))
[~,i_selected] = min(abs(phase - selected_phase));

P = length(phase);
const = NaN * ones(1,P);
amp = NaN * ones(1,P);
for j = 1:P
%  if exist('gamma') % DEBUG this is main Octave/Matlb version, just excluding "addi"
  % Evaluate a two column matrix representing the constant and oscillating term, for a least square fit by pseudoinverse (slash operator)
  A = [ones(N,1), cos(2*pi*time/period + phase(j)*pi/180)];
  v = A \ raw_signal; % solve linear equation system raw_signal = A v, v = [const; amp]
%  else % less accurate estimate for addi
%    v = [mean(raw_signal) ; mean((raw_signal-mean(raw_signal))./cos(2*pi*time/period + phase(j)*pi/180) )];
%  end
  const(j) = v(1);
  amp(j) = v(2);
end
all = struct('phase',phase, 'amp',amp, 'const',const);
[amp_max, i_max] = max(amp);
peak_phase = phase(i_max);
amp_rel_max = amp / amp_max;
%amp_rel_var = (amp-min(amp)) / (amp_max - min(amp)); % not used, as I want halfwidth = 360 deg if all amplitudes are almost equal
% and since opposite phase gives negative amplitude, it will be similar to 0.5*amp_rel_max (at least on noise free data)

% Find FWHM of the amp_rel_max array
%TODO: fwhm_phase = fwhmi(...

% Simpler estimate
is_above = amp_rel_max >= 0.5;
if sum(diff(is_above)>0) > 1
  error('Possibly more than one local maximum.')
end

i = 1:P;
i1 = find(~is_above(1:i_selected), 1, 'last')+1 % first is_above without gap from i_selected
i2 = find(~is_above(i_selected:P), 1, 'first')-2+i_selected % last is_above without gap from i_selected
% FIXME handle wrapping in case i1 or i2 is empty array



fullwidth_phase = (i2-i1) * search_step;
