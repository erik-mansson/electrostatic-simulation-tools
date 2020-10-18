% Model function to do a sum of multiple Gaussians.
%
% PARAMETERS
%  abscissa   The axis on which to evaluate the function(s).
%  areas      N-by-1 array with the area of each peak.
%             The area aparameter means that if an integral accountinf for bin width is made, then
%             sum(yield * abscissa_bin_width) = sum(areas).
%  widths     N-by-1 array with the FWHM of each peak.
%  centres    N-by-1 array with the central position of each peak.
% 
% EXAMPLES
% To fit area, width and centre of each peak:
% initial = [ones(length(masses),1), repmat(0.3,length(masses),1), masses'];
% model   = @(params,x) multiple_gaussians_model(x, params(:,1), params(:,2), params(:,3));
%
% To fit area and width of each peak:
% initial = [ones(length(masses),1), repmat(0.3,length(masses),1)];
% model   = @(params,x) multiple_gaussians_model(x, params(:,1), params(:,2), masses');
%
% To fit only area of each peak:
% initial = [ones(length(masses),1)];
% model   = @(params,x) multiple_gaussians_model(x, params(:,1), 0.3, masses');
% 
% Then fit using:
% result_info = struct();
% [values, resnorm, result_info.residuals, result_info.exitflag, various, ~, result_info.Jacobian] = ...
%     lsqcurvefit(model, initial, abscissa, signal, lower, upper, opt);
% [half_confint_68, sigma_t] = fitparam_std_estimate(values, result_info);
%
function [yield] = multiple_gaussians_model(abscissa, areas, widths, centres)

if size(abscissa,1) ~= 1
  error('abscissa should be a row vector.')
elseif size(areas,2) ~= 1
  error('areas should be a column vector.')
elseif size(widths,2) ~= 1
  error('widths should be a column vector.')
elseif size(centres,2) ~= 1
  error('centres should be a column vector.')
end
sigmas = widths / (2*sqrt(2*log(2))); % convert from FWHM to sigma-parameter

  
yield = sum(normpdf(abscissa, centres, sigmas) .* repmat(areas,1,length(abscissa)),1);
