%  yield = decay_model(param, t)
%
% Model for pump-probe yield given by the convolution of
% pulse duration(s) and the exponential decay of an intermediate state.
% But generally useable to calculate the convolution between
% a gaussian and an exponential decay (where the decay is multiplied by a
% heaviside function so that the decaying function is initially zero).
%
% In the process: initial state --XUV--> A --IR-----> B
%                         competing with A --decay--> C
% the B state is the one observed in the ion yield (or photoelectron) signal.
% A competing process is the decay of the A-state with lifetime_A to C
% which prevents the B-state from being reached via the IR transition.
%
% The remaining population in the A state and thus signals for
% probing A--IR-->B as well as fragments emitted automatically in the decay _to_ C
% will be proportional to
%   yield = convolution(convolution(I_XUV, P_A_not_decayed), I_IR)
% and the article rewrites this as
%   yield = convolution(convolution(I_XUV, I_IR), P_A_not_decayed)
%   yield = convolution(I_combined, P_A_not_decayed)
% where the IR laser and HHG XUV pulses have been assumed to be
% Gaussian in intensity:
%   I_XUV = I0_XUV * exp(- t^2 / (2 * std_XUV^2))
%   I_IR = I0_IR * exp(- t^2 / (2 * std_IR^2))
% giving the convolution of both intensities
%   I_combined = I0_combined = exp(- t^2 / (2 * std_combined^2))
%   std_combined = sqrt(std_XUV^2 + std_IR^2), "pump-probe resolution" (std.dev., not FWHM)
% The decay is exponential:
%   P_A_not_decayed = exp(-t / lifetime_A)
% To return this convolution result, let the lifetime_A be positive!
% The article then gives the final formula (12) for "yield"
% which is implemented here.
% Note: this model (Lacoursiere, Meyer, Nahon, Morin, Larzilliere; Nucl.Instr. 1994)
% seems to be somewhat simplistic or intended for long timescales
% (No attempts to handle fields, coherence and subcycle effects or nonlinearities.
%  Just incoherent "classical" probabilities & lifetimes are treated.)
%
% If the C state (after decay) is assumed to be stable (not decaying further on relevant time scale)
% then signal from probing C--IR-->D would be proportional to
%   yield = convolution(convolution(I_XUV, P_A_has_decayed), I_IR)
%         = convolution(convolution(I_XUV, 1 - P_A_not_decayed), I_IR)
% To return this convolution result, let the lifetime_A be negative!
% 
% This function is intended to be used by lsqcurvefit or similar optimizers.
%
% PARAMETERS
%   param  the model parameters, as a column vector:
%          [lifetime_A; std_combined; t_zero; background; factor]
%          NOTE: the global variables 
%             decay_model_fixed_std_combined
%             decay_model_fixed_t_zero
%             decay_model_fixed_background
%          will if nonempty replace the use of some parameter columns, so
%          that a smaller parameter matrix should be sent (indexing changes).
%          
%          To convert from FWHM: std_combined = FWHM_combined/(2*sqrt(2*log(2)))
%          
%          2014-03-05: Now also allowing negative lifetime_A.
%          * Positive lifetime means that the transient state is probed,
%            population (and signal) is decreasing with probe delay.
%            Called 'decay1' in mass_scan_analysis.m.
%            Shape: height*exp(-t/lifetime_A)
%          * Negative lifetime means that the state that decay leads to is probed,
%            population (and signal) is increasing with probe delay.
%            Called 'final' in mass_scan_analysis.m.
%            Shape: height*(1 - exp(-t/lifetime_A))
%          2016-10-25
%          * TODO (or another function): allow something like A--decay1->C--decay2->D
%            and IR-probing from C, so that thre is both a slow rise due to decay1
%            and a signal decay due to decay2 competing with IR-probing.
%
%   t      an array or matrix of physical times [fs]
%          for which the yield should be evaluated.
%          (t = scan_value should approximte the experiment).
%
% RETURNS
%   yield  an array or matrix of the same size as t, with the
%          calculated yield (probability of reaching the B-state)
%          at the times t. The unit is undefined, same as the scaling
%          factor = param(5);
% 
%          Assuming t contains a suffiently wide interval and that lifetime_A is positive:
%          * The result is normalized so that
%            sum(yield-background)*dt = factor*lifetime_A = param(5)*param(1).
%            which means the step height is equal to factor=param(5) in the limit
%            of std_combined=0 (before the convolution).
%          * The weighted average time is t_zero+lifetime_A = param(3)+param(1).
%          For negative lifetimes (not probing transient state), the integral of the result is not bounded.
%
function yield = decay_model(param, t)

% Model parameters
global decay_model_fixed_std_combined decay_model_fixed_t_zero decay_model_fixed_background

probing_before_decay = param(1) > 0; % was always true before 2014-03-05 (and actually lifetime_A was not guaranteed to be positive?)
lifetime_A   = abs(param(1)); % [fs] lifetime of the imtermediate state populated by the XUV pulse

% Some parameters can be held fixed by the use of global variables!
% This means the input param() array is shorter.
if length(decay_model_fixed_std_combined) == 0
  std_combined = param(2); % [fs] "pump-probe resolution" = sqrt(std_XUV^2 + std_IR^2)
  p = 3;
else % don't optimize the std_combined, use a fixed value
  std_combined = decay_model_fixed_std_combined;
  p = 2;
end
if length(decay_model_fixed_t_zero) == 0
  t_zero       = param(p); % [fs] the time when both pulses overlap, an experimental offset from zero.
  p = p + 1;
else % don't optimize the t_zero, use a fixed value
  t_zero       = decay_model_fixed_t_zero;
end
if length(decay_model_fixed_background) == 0
  background   = param(p); % same unit as the experimental signal that yield is compared to (e.g. [ions/laser shots] or [fraction of ion diagram = count(selected ion)/count(any ion)] or [count/s])
  p = p + 1;
else
  background   = decay_model_fixed_background;
end
factor       = param(p); % same unit as the experimental signal that yield is compared to (e.g. [ions/laser shots] or [fraction of ion diagram = count(selected ion)/count(any ion)] or [count/s])

if imag(t_zero)~=0
  error('Non-real t_zero: %s', num2str(t_zero));
end

% Get "true" delay, accounting for experimental offset
t = t - t_zero;

if std_combined <= 0
  % If no gaussian width, then don't do the convolution.
  if isinf(lifetime_A) && abs(lifetime_A) > 1E100
    % Infinite lifetime ==> just a simple step function
    yield = background + factor * heaviside(t);
%     % If there is a bin near t==0 (exactly), we can make normalization closer
%     % to unitary by setting exp(0) as 1/2 intead of 1. This agrees with the limit
%     % of the convolution (as defined below) when std_combined approaches zero.
    yield(t == 0) = background + factor / 2;
  else
    % Just return the exponential decay.
    if probing_before_decay % signal decreases with delay (for t>>std_combined)
      t(t < 0) = Inf; % to only get background before t_zero (as convolution with step function)
      yield = background + factor * exp(-t/lifetime_A);
      % If there is a bin at t==0 (exactly), we can make normalization closer
      % to unitary by setting exp(0) as 1/2 intead of 1. This agrees with the limit
      % of the convolution (as defined below) when std_combined approaches zero.
      yield(t == 0) = background + factor / 2;
    else % probing after decay: signal increases with delay (for t>>std_combined)
      yield = background + factor * (1 - exp(-t/lifetime_A));
      yield(t <= 0) = background;
    end
  end
  
else
  % The convolved funcion (normal case):
  % The division by two makes the result normalized so that
  % sum(y-background)*dt = factor (if probing_before_decay).
  yield = background + (factor/2) ...
    * exp(-t/lifetime_A + 0.5*(std_combined/lifetime_A)^2) ...
    .* erfc((-t/std_combined + std_combined/lifetime_A) / sqrt(2));
  % NOTE: the last factor can equivalently be written with erf instead of erfc, so that the t stands without minus sign:
  % .* (1 + erf((t/std_combined - std_combined/lifetime_A) / sqrt(2)));
  % This form perhaps makes it more clear that this factor is step-function like and 
  % runs from 0 at t --> -Inf to 2 at t --> +Inf (motivating the use of "factor/2" earlier in the expression).
  
  if lifetime_A == 0
    % 2013-04-25 adding also the possibility to get unconvolved Gaussian by
    % setting zero lifetime. In this case the normalization is different,
    % "factor" is the area, rather than area/lifetime_A.
    %yield = background + ( factor / sqrt( 2*pi*std_combined^2) ) ...
    %  * exp(- t.^2 / (2 * std_combined^2)); % equivalent to:
    yield = background + factor * normpdf(t, 0, std_combined);
    
    % When using in mass_scan_anlysis or any kind of fitting, this behaviour in the limit of zero lifetime makes sense
    % for fitting the model directly,
    % but NOT when the curve expression involves the difference between two decay_model-calls (like for 'decay2', 'decay2t' and 'decay2a' in mass_scan_analysis.m)
    warning('decay_model:zero_lifetime', 'Returning a pure Gaussian as special case when zero lifetime_A given.');
    
  elseif (-min(t)/lifetime_A + 0.5*(std_combined/lifetime_A)^2) > 698
    % Discovered 2012-01-05 that when very wide compared to lifetime,
    % i.e. similar to a symmetric gaussian, the numeric evaulation breaks
    % down because exp(0.5*(std_combined/lifetime_A)^2) becomes Inf (but exp(709) is 8.2E307).
    % Since the decay is so short in this case, it should be OK to
    % approximate the function as the unconvolved gaussian.
    % To follow the same normalization as in the convolved case,
    % a multiplication by lifetime_A is done.
    % To get the "center of intensity" at t_zero+lifetime_A,
    % the Gaussian is shifted by lifetime_A.
    y = background + ( factor * lifetime_A / sqrt( 2*pi*std_combined^2) ) ...
      * exp(- (t-lifetime_A).^2 / (2 * std_combined^2));
    
    % Replace these points:
    exponent = -t/lifetime_A + 0.5*(std_combined/lifetime_A)^2;
    which = isnan(yield) | exponent>699;
    yield(which) = y(which);
    % Interpolate these points (to soften transition):
    which = exponent > 696;
    weight = min(1, (exponent(which) - 696)/(699 - 696));
    weight = weight.^2 ./ ((1-weight).^2 + weight.^2); % make the derivative of weight and yield continuous, to make the interpolation less notable
    yield(which) = yield(which).*(1-weight) + weight.*y(which);
  end

  if ~probing_before_decay % probing after decay: signal increases with delay (for t>>std_combined)
    % Adding this ad-hoc after the "probing before decay"-signal has been computed.
    decay_convolved = yield - background;
    gaussian_CDF = factor * normcdf(t, 0, std_combined);
    % Not sure how it was derived, but in the limit of std_combined-->0
    % this expression approaches the desired:
    %   yield ={ background + height*(1-exp(-(t-t0)/lifetime_A)), t>t0
    %          { background                                     , t<t0
    yield = background + gaussian_CDF - decay_convolved;
    
    % NOTE: the result for lifetime_A = -Inf becomes flat (zero height)
    % but it makes sense, since +Inf gives a non-decaying step so no 
    % population would reach a final state.
  end
  
end

% Note: according to the article's equation (13)
% excitation_A_rate * (peak_IR_flux * cross_section_A_to_B) * detection_efficiency =
%   = C * sqrt(2) std_combined  / (std_XUV * std_IR * sqrt(pi))

% TODO: test cftool