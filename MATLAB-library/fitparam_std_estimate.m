% Estimate standard deviation in parameters found by an optimizer
% such as fit, lsqnonlin or lsqcurvefit.
% [fitted,~,info] = fit(...)
% or
% [fitted,~,info.residuals,~,~,~,info.Jacobian] = ...
% PARAMETER
%  fitted     Result of fit as a cfit object or an array of parameter values.
%             If fitted is a matrix (possible with lsqcurvefit),
%             then it will be converted to array form as fitted(:).
%  info       Structure with the fields 'residuals' and 'Jacobian'.
%             If residuals is a matrix (possible with lsqcurvefit),
%             then it will be converted to array form as residuals(:).
%  may_use_curvefit_confint (Default: true) Call confint() if a cfit-object
%             is given, assuming that the algorithm there uses all
%             available info (Rinv there is probably more like the pseudoinverse 
%             of the "basis matrix" in my linear least squares formulation,
%             and Jacobian should be available too if useful.)
%             For DI2nothigh on u2e the confint gives a few percent larger
%             value so I prefer to use that (may_use_curvefit_confint=true.
% RETURNS
%  half_confint_68    Estimated halfwidth of 68% confidence intervals (same confidence level as standard deviation in case of normal distribution).
%                     This was previously called "sigma" although it for Student's t-distribution is always somewhat larger (depending on degrees_of_freedom).
%                     Using it like sigma will still mean that differences relative to this are tested at the 68% confidence level.
%                     To compute a correct p-value, any signal or difference, direct division by half_confint_68 is not correct, need to use the degrees_of_freedom
%                     to convert it to the sigma_of_t_distribution.
%  sigma_t            The actual "sigma" value of the Student's t-distribution.
%                     A two-sided confidence interval for some p_value is obtained as:
%                     fitted_value +- sigma_t * tinv(1 - p_value/2, degrees_of_freedom).
%                     For any relative difference between fitted parameters
%                       t = (value1-value2)/sqrt(sigma_t1^2 + sigma_t2^2)
%                       is t-distributed with
%                       p_vaule_twosided = (1-tcdf(t, nu))*2 % =tcdf(-t, nu)+1-tcdf(t, nu)
%                       p_value_onesided =  1-tcdf(t, nu) 
%                     Since many programs are showing the half_confint_68 instead of sigma_t, one needs to do a rescaling to get
%                     the proper t and p-values for a test (or confidence intervals at other confidence levels) :
%                       t_incorrect = (value1-value2)/sqrt(half_confint_68_1^2 + half_confint_68_2^2) 
%                       t_incorrect = (value1-value2)/sqrt(sigma_t_1^2*tinv(1-alpha32/2, degrees_of_freedom1)^2 + sigma_t_2^2*tinv(1-alpha32/2, degrees_of_freedom2)^2) 
%                       In case both parameters are from the same fit (same degrees_of_fredom) it can be easily rescaled afterwards:
%                       t_incorrect = (value1-value2) / sqrt(sigma_t_1^2 + sigma_t_2^2) / tinv(1-alpha32/2, degrees_of_freedom1)
%                       ==>
%                       t = t_incorrect * tinv(1-alpha32/2, degrees_of_freedom1)
%                     then the usual equations give p-values
%                       p_vaule_twosided = (1-tcdf(t, nu))*2 % =tcdf(-t, nu)+1-tcdf(t, nu)
%                       p_value_onesided =  1-tcdf(t, nu) 
%
%  half_confint_95    Estimated halfwidth of 95% confidence intervals.
%  degrees_of_freedom = number of datapoints - number of fitted parameters
%                     needed to compute other p-values or confidence intervals using Student's t-distribution.
%
% SEE ALSO
%  confint which may be used after fit, with alpha as second argument.
%
function [half_confint_68, sigma_t, half_confint_95, degrees_of_freedom] = fitparam_std_estimate(fitted, info, may_use_curvefit_confint)
if nargin < 3
  may_use_curvefit_confint = true;
end

alpha32 = (1-normcdf(1))*2; % p_value=approx=0.32 (twosided) to give approximately 68% confidence interval, assuming normal distribution

original_fitted_shape = size(fitted);
done = false;
if strcmp(class(fitted), 'cfit')
  % When fit() was used
  beta = coeffvalues(fitted);
  res = info.residuals;
  p = length(beta);
  
  degrees_of_freedom = numel(res) - numel(beta); % degrees of freedom
  
  if may_use_curvefit_confint
    % 2013-07-04: assuming the confint(cfit) can do best possible estimate.
    % It uses cftinv insted of tinv for t-inverse, and names are different
    % but the algorithm (Rinv of something like a basis matrix, seems
    % related to my use of pseudoinverse for the linear equation system
    % variant.)
    try
      half_confint_68 = diff(confint(fitted,1-alpha32))/2;
      half_confint_95 = diff(confint(fitted,1-0.05))/2;
      
      %very old sigma_from95 = diff(confint(fitted,1-alpha32))/2 / norminv(1-0.05/2);
      % Since a t-distribution is used in confint, the underlying sigma can be obtained by dividing out the scaling factor used for the chosen p-value level:
      %sigma_from95 = half_confint_95 / tinv(1-0.05/2, numel(res)-p); % equivalent used before 2019-04-07 when starting to use the variable name degrees_of_freedom
      sigma_from68 = half_confint_68 / tinv(1-alpha32/2, degrees_of_freedom);
      sigma_from95 = half_confint_95 / tinv(1-0.05/2, degrees_of_freedom);
      
      done = true;
    catch e
      warning('fitparam_std_estimate:nlparci', 'Could not determine confidence intervals using confint(), will try nlparci...');
      e
    end
  end

  % previously:
  J = info.Jacobian;
elseif isnumeric(fitted)
  % Assume a structure constructed with info arguments from lsqnonlin
  if numel(fitted) > length(fitted)
    % If fitted is a matrix (possible with lsqcurvefit),
    % then it will be converted to array form as fitted(:).
    fitted = fitted(:);
  end
  beta = fitted;
  p = numel(beta); % number of parameters fitted

  if isfield(info,'residuals')
    res = info.residuals;
  elseif isfield(info,'residual')
    res = info.residual;
  else
    error('No field named ''residuals'' found in fit info structure.')
  end
  if numel(res) > length(res)
    % If residuals is a matrix (possible with lsqcurvefit),
    % then it will be converted to array form as residuals(:).
    res = res(:);
  end
  degrees_of_freedom = numel(res) - numel(beta); % degrees of freedom
  if isfield(info,'Jacobian')
    J = info.Jacobian;
  elseif isfield(info,'jacobian')
    J = info.jacobian;
  else
    error('No field named ''Jacobian'' found in fit info structure.')
  end
  if size(J,1) ~= numel(res) || size(J,2) ~= p
    error('Unexpected size relation between dataset residuals, parameter array and Jacobian matrix.'); % NOTE: this check is new so there be some case where it falsely raises the error...
  end
end
if issparse(J)
  r = sprank(J); % probably an upper bound for rank
  if r >= p % rank may be OK, check using full calculation too
    J = full(J);
    r = rank(J);
  end
else
  r = rank(J);
end
if r < p
  warning('fitparam_std_estimate:rank', 'Jacobian does not have full rank, errorbars may be invalid.')
end
if numel(res) <= p
  warning('fitparam_std_estimate:underdet', 'Too many parameters for the dataset, errorbars are invalid.')
end

if ~done
  
% (Original from 2014):
% (1-alpha)*100% confidence intervals are [mu-delta to mu+delta] where
% delta = z*sigma and z is the conversion factor. Thus sigma = diff(conf)/2/z.
% For a normal distribution: alpha=0.05 gives z=1.96, alpha=alpha32 gives z=1.
% The nlparci uses Student's t-distribution rather than normal
% distribution, to more correctly account for the number of estimated
% parameters. (But still I think it is a gross approximation because the
% motivations I have seen for t-distribution (p=1 case) discuss that a
% parameter is estimated using n measurements, while here we are handling
% the case where an optimization of the parameter used n points on a
% function curve. For histogram data, binning it differently will yield
% different n.
% For n>>p and noise in each bin normally distributed with standard deviation
% sigma_exp, where the model IS ABLE to perfectly fit it is perhaps still
% reasonable to think the noise dominates the residual, from which
% sigma_exp would be estimated. However, if the model never fits very well
% or if the noise is Poissonian there will be modeling error.
% What I don't yet understand is how the z that seems to be derived when discussing
% direct noise sigma_exp in the measurement, is by nlparci assumed to be
% relevant also for uncertainty in the parameter estimates.
%
% 2019:
% Well, it makes sense that (degrees of freedom = ) number of datapoints minus number of estimated parameters as input
% should influence the uncertainty (motivating Student-t instead of Normal distribution, converging for large number of points and few parameters).
%
% The z-value is shared by all parameters in nlparci, and as it is able to
% discover when a few parameters are redundant while others are
% well-determined it is obviouly the estimation of sigma* (for each
% parameter) that is the important part.
% Algorithm seems to agree with http://www.gnu.org/software/gsl/manual/html_node/Computing-the-covariance-matrix-of-best-fit-parameters.html
% although completely different notation.

  half_confint_68 = NaN(size(beta));
  sigma_from95    = NaN(size(beta));
  try
    half_confint_68 = diff(nlparci(beta,res,'jacobian',J,'alpha',alpha32),[],2) / 2;
    half_confint_95 = diff(nlparci(beta,res,'jacobian',J,'alpha',0.05),[],2) / 2;
    
    %sigma_from95    = half_confint_95 / norminv(1-0.05/2); % NOTE: was in used until 2019-04-07 (although t-distribution was already used in the case that handles cfit result objects!)
    % Since a t-distribution is used in confint, the underlying sigma can be obtained by dividing out the scaling factor used for the chosen p-value level:
    sigma_from68 = half_confint_68 / tinv(1-alpha32/2, degrees_of_freedom);
    sigma_from95 = half_confint_95 / tinv(1-0.05/2, degrees_of_freedom);

  catch e
    warning('fitparam_std_estimate:nlparci', 'Could not determine confidence intervals, possibly too large Jacobian');
    e
  end

end



% % Development tests:
% % nlinfit, nlparci
% c95 = nlparci(coeffvalues(fitted),info.residuals,'jacobian',info.Jacobian,'alpha',0.05);
% c68 = nlparci(coeffvalues(fitted),info.residuals,'jacobian',info.Jacobian,'alpha',(1-normcdf(1))*2);
%  % z = norminv(1-0.05/2) == 1.9600 <==> alpha=(1-normcdf(1.96))*2 == 0.0500 (p=0.9500)
%  % one half_confint_68 (z=1) corresponds to alpha=(1-normcdf(1))*2=0.3173 (p=0.6827)
% %s95f = diff(confint(fitted)',[],2) / norminv(1-0.05/2);
% % To almost four digits s95 and s68 agree (ratio 1.0002)
% s95 = diff(c95,[],2)/2 / norminv(1-0.05/2)
% s68 = diff(c68,[],2)/2

% Old, used until 2019-04-07
a = (sigma_from95 + half_confint_68)/2;
a = a + eps(a) * 1E2;
if max(abs(sigma_from95-half_confint_68) ./ a) > 0.04
  warning('fitparam_std_estimate:differ', 'Sigma computed from 68%% and 95%% confidence intervals differ by more than 4%%. (n=%d,p=%d)\n    s68=%s\n    s95=%s', length(res), p, mat2str(half_confint_68',4), mat2str(sigma_from95',4));
end

% 2019-04-19 after fixing to use Student's t-distribution and compare with real underlying sigma (estimate) rather than 68% confidence interval (estimate). Can use much stricter limit then.
a = (sigma_from95 + sigma_from68)/2;
sigma_t = a; % use the average in case of any numeric noise
a = a + eps(a) * 1E2;
if max(abs(sigma_from95-sigma_from68) ./ a) > 0.0001
  warning('fitparam_std_estimate:diagree', 'Sigma computed from 68%% and 95%% confidence intervals differ by more than 0.01%%. (n=%d,p=%d)\n    s68=%s\n    s95=%s', length(res), p, mat2str(half_confint_68',4), mat2str(sigma_from95',4));
end

% 2019:
% The Student's t-test or confidence interval using tcdf() tinv()
% would use t = value/sigma, or for a two-distribution test: t = (value1-value2)/sqrt(sigma1^2 + sigma2^2)
% and the degrees of freedom nu = numel(data) - numel(parameters_fitted)
%   https://en.wikipedia.org/wiki/Student's_t-distribution#In_frequentist_statistical_inference
%   https://en.wikipedia.org/wiki/Degrees_of_freedom_(statistics)
% nu = numel(output.signal_datapoints) - numel(output.internal.values); % degrees of freedom
% p_vaule_twosided = (1-tcdf(t, nu))*2 % =tcdf(-t, nu)+1-tcdf(t, nu)
% p_value_onesided =  1-tcdf(t, nu) 
% E.g.
% 1-68% =approx= alpha32 = (1-tcdf(t_alpha32, nu))*2
%  <==> tcdf(t_alpha32, nu) = 1 - alpha32 / 2 
%  <==> t_alpha32 = tinv(1 - alpha32/2, nu)  ; % e.g. 1.0906 for nu=6 as an example, i.e. a bit more than the 1 that Normal distribution would give (norminv(1 - alpha32/2, 0, 1))
% t_alpha5 = tinv(1 - 0.05/2, nu)  ; % e.g 2.4469 for nu=6 as an example, i.e. a lot more than the 1.96 that Normal distribution would give (norminv(1 - 0.05/2, 0, 1))

% To answer whether any conversion is needed of the sigma (s68) returned by this function, in order to compute other confidence intervals.
% As a first test of which distribution (and nu if Student's t) was used:
%  confint95/confint68 =  1.9600 / 1 for normal distribution
%  confint95/confint68 = tinv(1 - 0.05/2, nu) / tinv(1 - alpha32/2, nu)  % = 2.2437 for nu = 6
%  confint95/confint68 = tinv(1 - 0.05/2, nu) / tinv(1 - alpha32/2, nu)  % = 1.9836 for nu = 60
%  confint95/confint68 = tinv(1 - 0.05/2, nu) / tinv(1 - alpha32/2, nu)  % = 1.9717 for nu = 120, i.e. now quite close to normal distribution
% Example:
%   nu = numel(output.signal_datapoints) - numel(output.internal.values)    % = 92
%   output.bg.conf95halfrange ./ output.bg.std  % =  1.9753    1.9753    1.9753    1.9753    1.9753    1.9753
%   % the same ratio by tau_B, tau_A (for masses where fitted), and the common t0 -- expected since same nu for all parameters in a fit.
%   tinv(1 - 0.05/2, nu) / tinv(1 - alpha32/2, nu) % == 1.9753
%   I.e. perfectly matched ratio, confirms that Student's t distribution was used and that the _std is not actually like "sigma" but a few percent larger (dependent on nu)
%
% To compute p-value for a Student's t-test, one needs the test statistic
%   t = value/real_sigma   or for a two-distribution test: t = (value1-value2)/sqrt(sigma1^2 + sigma2^2).
% where real_sigma can be obtained as
%   sigma_from_t95 = conf95halfrange/tinv(1 - 0.05/2, nu)
%   sigma_from_t68 = conf68halfrange/tinv(1 - alpha32/2, nu)
% 
% sigma 

% s68 = diff(c68,[],2)/2 = 
% c68 = nlparci(coeffvalues(fitted),info.residuals,'jacobian',info.Jacobian,'alpha',(1-normcdf(1))*2);

if any(size(fitted) ~= original_fitted_shape)
  % If fitted is a matrix (possible with lsqcurvefit),
  % then it will be converted to array form as fitted(:) during processing,
  % but we can return uncertainty values in the original shape.
  half_confint_68 = reshape(half_confint_68, original_fitted_shape);
  sigma_t         = reshape(sigma_t,         original_fitted_shape);
  half_confint_95 = reshape(half_confint_95, original_fitted_shape);
end

  
