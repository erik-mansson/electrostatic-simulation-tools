% iif
% Inline if-expression.
% If the condition evaluates to true (not false or 0)
% true_value is returned. Otherwise false_value is returned.
%
% Any parameters may be matrix or scalar. All matrix parameters
% must be of the same size.
%
% PARAMETERS
%  condition    decides what to return [boolean]
%  true_value   what to return if all(not(condition))
%  false_value  what to return if all(condition)
% RETURN
%  If condition is scalar then true_value or false_value without change.
%  If condition is a matrix then the result is a matrix of that size.
%  In case of scalar values they are repeated to fill the required size.
function result = iif(condition, true_value, false_value)
if length(condition) == 1
  % Scalar condition, possibly non-scalar values:
  if condition
    result = true_value;
  else
    result = false_value;
  end

elseif length(true_value) == 1
  if length(false_value) == 1
    % Matrix condition. Scalar values.
    result = ones(size(condition)) * false_value;
    result(condition ~= false) = true_value;
  elseif all(size(false_value) == size(condition))
    % Matrix condition and false_value. Scalar true_value.
    result = false_value;
    result(condition ~= false) = true_value;
  else
    error('Unsupported parameter size combination: %s, %s, %s.', ...
      mat2str(size(condition)), mat2str(size(true_value)), mat2str(size(false_value)))
  end

elseif all(size(true_value) == size(condition))
  if length(false_value) == 1
    % Matrix condition and true_value. Scalar false_value.
    result = true_value;
    result(not(condition)) = false_value;
  elseif all(size(false_value) == size(condition))
    % Matrix condition, true_value and false_value (all of same size).
    result = true_value;
    false_indices = not(condition);
    result(false_indices) = false_value(false_indices);
  else
    error('Unsupported parameter size combination: %s, %s, %s.', ...
      mat2str(size(condition)), mat2str(size(true_value)), mat2str(size(false_value)))
  end

else
  error('Unsupported parameter size combination: %s, %s, %s.', ...
    mat2str(size(condition)), mat2str(size(true_value)), mat2str(size(false_value)))
end