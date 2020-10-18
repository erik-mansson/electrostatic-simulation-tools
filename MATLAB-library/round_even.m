% An implementation of unbiased rounding to integers, with ties towards even.
% 
% This is equivalent to the convergent() function, but much faster (about 55 times faster on 1280x1024 matrices of doubles).
% The difference occurs for -0.5, 0.5, 2.5, etc.
% which round() would return as -1, 1, 3 (away from zero)
% but unbiased ties-towards even rounding gives 0, 0, 2.
%
% SEE ALSO
%  convergent()  Equivalent but much slower (1280x1024 matrix of doubles in 3.9 s, instead of 0.07 s)
function rounded = round_even(values)

fraction = rem(values, 2); % reminder after division by 2: 0 <= rem < 2 for positive inputs, similar with negative sign for negative inputs
% In most cases round() does what we want, but we need to consider the cases where fraction is +-0.5 or +-1.5.
%   [values; round(values); convergent(values); fraction]
% round() gives the correct result for 1.5 and 3.5, since the higher 2 and 4 are even.
% round() fails for 0.5 and 2.5 (giving 1 and 3) which are 0 and 2 by convergent() round-to-even.
% Thus the failures for positive inputs are precisly those where rem(values, 2) == 0.5.
% For negative values the same relation applies with a minus sign in both values and fraction,
% meaning that abs(fraction) == 0.5 identifies the values that incorrectly would be rounded
% to larger absolute values. For positive values we can subtract somehting small to make them
% closer to the correct integer below. For negative values we want to add something small,
% which is equivalent to subtracting something small negative, i.e subtracting the fraction always works,
% for the matrix elements where abs(fraction) == 0.5.

% Zero all the fractions where no corection is needed.
fraction(abs(fraction) ~= 0.5) = 0;

% Apply correction and call round()
rounded = round(values - fraction);

% Compare results:
% values = -4.5:4.5; 
% values = -5:0.25:5;
% [values; round(values); convergent(values); rounded]