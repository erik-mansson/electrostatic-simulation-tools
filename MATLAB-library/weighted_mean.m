% Helper function for calculating weighted & unweighted averages, 
% with various uncertainty estimates for errors of the mean.
%
% The weights will be 1 ./ std_estimates.^2.
% Different attempts to quantify uncertainty in the average are also returned.
% 
% PARAMETERS
%  values         N-element array or (N-by-C)-matrix of values to average.
%  std_estimates  N-element array or (N-by-C)-matrix of estimated standard deviations for each value.
%                 NOTE: if any of the std_estimates is very close to zero, its weight will become dominat.
%                 You could consider using max(std_estimates, std(values)) to put a lower bound
%                 to the uncertainty estimates from a fitter, using the spread of the values.
%                 This is similar to the combination used in stderr_avg_combined, but by applying it to the input
%                 it also concerns the average that will be computed.
%
% RETURNS, all as (1-by-C) matrices/arrays
%  weighted_average                   Weighted average.
%
%  std_sample_unbiased                Unbiased estimated of sample variance, after accounting for the weights.
%
%  stderr_avg_combined                Intended to indicate uncertainty in weighted_average,
%                                     stderr_avg_combined = max([
%                                                 sqrt( (stderr_avg_by_propagation^2 + stderr_avg_by_sample_div_by_sqrtN^2)/2 )
%                                                 stderr_avg_by_propagation
%                                     ]);
%                                     The reason for using stderr_avg_by_propagation as a lower bound is that
%                                     stderr_avg_by_sample_div_by_sqrtN is not well-motivated, but still seems interesting to use to take
%                                     some noitice in case the sample is very spread although each point has small uncertainty
%                                     (which would question whether the values are really from a distribution with a common expectation value).
%
%  stderr_avg_by_propagation          Uncertainty in weighted_average given only by std_estimates (not actual variation between values).
%                                     This should be a standard error of the computed weighted_average, if one completely trusts that
%                                     there is no actual variation between each data point's expectation value, just random deviations
%                                     due to noise matching the estimated uncertaintites.
%
%  stderr_avg_by_sample_div_by_sqrtN  Rescaled version of std_sample_unbiased, dividing by sqrt(number of nonzero weights) 
%                                     although this is not strictly telling anything about sample mean uncertainty.
%
%  unweighted_average                 Unweighted average = mean(values).
%
%  unweighted_std_sample                Unweighted sample standard deviation = std(values).
%
%  unweighted_stderr_avg_by_propagation  Uncertainty in weighted_average given only by std_estimates (not actual variation between values).
%
% Author: Erik Månsson, erik.maansson@desy.de
function [weighted_average, std_sample_unbiased, ...
          stderr_avg_combined, stderr_avg_by_propagation, stderr_avg_by_sample_div_by_sqrtN, ...
          unweighted_average, unweighted_std_sample, unweighted_stderr_avg_by_propagation] ...
        = weighted_mean(values, std_estimates)


if min(size(std_estimates)) ~= 1
  error('The std_estimates must be a row or column array.');
end

if size(values,1) ~= length(std_estimates)
  if size(values,1) == 1 && size(values,2) == length(std_estimates)
    % OK to transpose one-dimensional arrays
    std_estimates = std_estimates(:); % make it a column array
    values        = values(:);        % make it a column array
  end
end
if ~all(size(std_estimates) == size(values))
  error('Mismatching size of values and weights.');
end

% According to https://en.wikipedia.org/wiki/Weighted_arithmetic_mean, using 1/std^2 as weight gives 
% "the maximum likelihood estimator of the mean of the probability distributions under the assumption
% that they are independent and normally distributed with the same mean."
weights = std_estimates .^ (-2);

invalid = ~isfinite(std_estimates) | std_estimates < 0;
if any(invalid)
  warning('weighted_mean:exclusion', '%d values excluded from weighted_mean, at indices %s.', sum(sum(invalid)), mat2str(find(invalid)));
end
weights(invalid) = 0; % 2018-11-12: now zeroing the weights of invalid points, making sure they don't affect the weighted mean


unweighted_average    = mean(values);
unweighted_std_sample = std(values);

% Ensure sum of weights is normalized to 1 within each column
weights = weights ./ repmat(sum(weights), size(weights,1), 1);


weighted_average                   =       sum(weights    .* values          , 1);
% The standard error of the computed weighted_average, can be computed using only "error propagation"
%   Var(sum(w_i * value_i)) = sum(Var(w_i * value_i)) =approx= sum(w_i^2 * std_estimate_i^2).
% This should be a standard error of the computed weighted_average, if one completely trusts that
% there is no actual variation between each data point's expectation value, just random deviations
% due to noise matching the estimated uncertaintites.
stderr_avg_by_propagation          = sqrt( sum(weights.^2 .* std_estimates.^2, 1) ); % This is my own approach, no reference except general https://en.wikipedia.org/wiki/Propagation_of_uncertainty

% Error propagation can be applied also to the unweighted mean's standard error
unweighted_stderr_avg_by_propagation  = sqrt( sum((1/length(values)).^2 .* std_estimates.^2, 1) ); % This is my own approach, no reference except general https://en.wikipedia.org/wiki/Propagation_of_uncertainty


% Like sample standard deviation, but defined in relation to weighted_average.
% https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance, assuming "Reliability weights" for unbiasing the estimate
%std_sample_biased                 = sqrt( sum(weights    .* (repmat(weighted_average,size(values,1),1) - values).^2 ));
std_sample_unbiased                = sqrt( sum(weights    .* (repmat(weighted_average,size(values,1),1) - values).^2) / (1-sum(weights.^2)/1) );

% If the same bias-correction is applied to my stderr_avg_by_propagation
%weighted_avg_err_by_propagation_u = sqrt( sum(weights.^2 .*  std_estimates.^2                                 ) / (1-sum(weights.^2)/1) );

% WARNING: this is not well-motivated, dividing a sample standard deviation by sqrt(number of values). If the average was not weighted, it would however be an estimate of error of the mean assuming uniform noise
stderr_avg_by_sample_div_by_sqrtN  = std_sample_unbiased ./ sqrt(sum(weights>0)); 

stderr_avg_combined = sqrt( (stderr_avg_by_propagation^2 + stderr_avg_by_sample_div_by_sqrtN^2)/2 ); % 2-norm-average
% Note that stderr_avg_combined can be smaller than underlying uncertainty estimates / sqrt(N), down to a factor 1/sqrt(2) in case the sample spread would be zero.
% Since only stderr_avg_by_propagation is well-motivated, it seems more proper to put it as lower bound:
stderr_avg_combined = max([stderr_avg_by_propagation stderr_avg_combined]);

