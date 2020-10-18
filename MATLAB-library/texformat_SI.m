% str = texformat_SI(value, significant_digits)
% Turn a number into a string with SI unit prefix,
% for use where TeX code is allowd, e.g. in figure legends and titles.
%
% PARAMETERS
%  value               The number to format into a string
%  significant_digits  The number of significant digits (Default: 3)
%                      Trailing zeroes are skipped, which may make fewer
%                      digits be returned. To keep trailing zeroes
%                      give a negative value for significant_digits.
%  use_greek_mu        Use TeX '\mu' rather than 'u' for micro (Default: true)
%
% RETURNS
%  string              Output string
% EXAMPLE 
%  texformat_SI(9.8765E-7)      % --> 0.99 \mu
%  texformat_SI(9.8765E-7, 4)   % --> 987.6 n
%  texformat_SI(9.8765E-7, 2)   % --> 0.99 \mu
%  texformat_SI(9.8765E7, 2)    % --> 9.9 M
%  texformat_SI(9.8765E7, 1)    % --> 0.1 G
%  texformat_SI(0.1,  2)       % --> 0.1 m
%  texformat_SI(0.1, -2)       % --> 0.10 m
%  texformat_SI(0.1,  4)       % --> 100 m
%  texformat_SI(0.1, -4)       % --> 100.0 m
%
% Created: 2012-02-27
% See also texformat_g10
% Author: Erik Månsson, Lund University, erik.mansson@sljus.lu.se
%
function str = texformat_SI(value, significant_digits, use_greek_mu)
keep_trailing_zeros = false;
if nargin < 2
%   if 10 ^ (3*mod(log10(abs(value))/3,1)) >= 500
%     % If the end would be >= 0.5 in some prefixed unit, enough to show two
%     % digits so that it becomes "0.5 G" rather than "500 M".
%     significant_digits = 2;
%   else
    significant_digits = 3;  
%   end
else
  if significant_digits < 0
    keep_trailing_zeros = true;
    significant_digits = max(1, round(-significant_digits));
  else
    significant_digits = max(1, round(significant_digits));
  end
end
if nargin < 3
  use_greek_mu = true;
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
        str = [str iif(j>1,', ','') texformat_SI(value(i,j), significant_digits)];
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
  v = -value;
  str = '-';
else
  v = value;
  str = '';
end
p = floor(log10(v));
s = floor((p+3 - min(3,significant_digits)) / 3);

if abs(s) < 1 || v == 0
  % No prefix
  if (significant_digits == 1 && v > 1) || keep_trailing_zeros
    % This write the desired number of digits, and always works,
    % but is is normally not preferred since we want to use "g" to skip
    % trailing zeroes. However, %.1g fails on 9.999 with significant_digits==1
    % which becomes 1e+001 rather than just 10. So when significant_digits==1
    % this method is used.
    str = [str sprintf(['%.' num2str(significant_digits-p-1) 'f'], v)];
  else
    str = [str sprintf(['%.' num2str(significant_digits) 'g'], v)];
  end
  % Appending space also when no prefix, for consistent appearance when unit symbol follows the output
  str = [str ' '];
  
elseif abs(s) > 8
  % Too large or small number, fall back to 10^p form
  str = texformat_g10(value, significant_digits);
  
else
  % Use some prefix
  if s > 0
    prefix = 'kMGTPEZY';
    prefix = prefix(s);
  else
    prefix = 'munpfazy';
    prefix = prefix(-s);
    if use_greek_mu && strcmp(prefix, 'u')
      prefix = '\mu';
    end
  end
  
  v = v/10^(3*s);
  if (significant_digits == 1 && v > 1) || keep_trailing_zeros
    % This write the desired number of digits, and always works,
    % but is is normally not preferred since we want to use "g" to skip
    % trailing zeroes. However, %.1g fails on 9.999 with significant_digits==1
    % which becomes 1e+001 rather than just 10. So when significant_digits==1
    % this method is used.
    str = [str sprintf(['%.' num2str(significant_digits-floor(log10(v))-1) 'f '], v) prefix];
  else
    str = [str sprintf(['%.' num2str(significant_digits) 'g '], v) prefix];
  end
end
