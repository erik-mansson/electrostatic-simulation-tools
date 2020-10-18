% str = texformat_g10(value, significant_digits)
% Turn a number into a string for use where TeX code is allowd,
% e.g. in figure legends and titles.
%
% PARAMETERS
%  value               The number to format into a string
%  significant_digits  The number of significant digits (Default: 2)
%  exponent_format     Alows customization of the exponent format.
%                     (Default: '\\cdot10^{%d}')
%                      Note that \ for TeX commands need double-escaping.
% RETURNS
%  string              Output string
% EXAMPLE 
%  texformat_g10(9.8765E-7)                      % --> 9.9\cdot10^{-7}
%  texformat_g10(9.8765E-7, 4)                   % --> 9.877\cdot10^{-7}
%  texformat_g10(9.8765E-7, 2, 'E%+d')           % --> 9.9E-7
%  texformat_g10(9.8765E-7, 2, 'E%+d')           % --> 9.9E-7
%  texformat_g10(9.8765E-7, 2, '\\times10^{%d}') % --> 9.9\times10^{-7}
%
% Created: <2009-01-06.
% Revised 2009-08-10 to handle negative values.
% Revised 2011-12-16 to allow customized exponent format.
% Author: Erik Månsson, Lund University, erik.mansson@sljus.lu.se
%
function str = texformat_g10(value, significant_digits, exponent_format)
if nargin < 2
  significant_digits = 2;
else
  significant_digits = round(significant_digits);
end
if nargin < 3
  exponent_format = '\\cdot10^{%d}';
end

if length(value) > 1
  if length(value) > 8
    error('Only scalar value supported, and arrays of at most 8 elements');
  else
    % Print array/matrix of values
    str = '';
    if size(value,2) == 1 && size(value,1) > 1
      % Turn column vector into row vector to use comma rather than
      % semicolon when one-dimensional array.
      value = transpose(value);
    end
    for i = 1:size(value,1)
      if(i>1)
        str = [str '; '];
      end
      for j = 1:size(value,2)
        str = [str iif(j>1,', ','') texformat_g10(value(i,j), significant_digits, exponent_format)];
      end
    end
  end
  return;
elseif length(value) == 0
  str = '[]';
  return;
elseif ~isfinite(value)
  str = sprintf('%g', value);
  return;
end

if value < 0
  value = -value;
  str = '-';
else
  str = '';
end
p = floor(log10(value));
needs_exponent = (p < -1 || p >= 3);

if ~needs_exponent && value >= 1000*(1 - 10^(-significant_digits)/2)
  % Value will be rounded up and then cross the p>=3 criterium for using
  % exponent notation (1E3 instead of 1000).
  needs_exponent = true;
  p = p+1;
end

if needs_exponent && value ~= 0
  str = [str sprintf(['%.' num2str(significant_digits) 'g' exponent_format], ...
                value/10^p, p)];
else
  str = [str sprintf(['%.' num2str(significant_digits) 'g'], value)];
end
