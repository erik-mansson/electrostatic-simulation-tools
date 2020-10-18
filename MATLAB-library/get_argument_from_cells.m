% value = get_argument_from_cells(options, name) % as 'bool'
% value = get_argument_from_cells(options, name, 'bool')
% value = get_argument_from_cells(options, name, 'val', default_if_not_named)
% value = get_argument_from_cells(options, name, 'opt', default_if_not_named, default_if_no_following_value)
%
% Helper function to use in other functions where a cell array of named
% arguments is accepted as input.
%
% PARAMETERS
%  options     Cell array of strings, possibly with non-strings after their
%              argument name string.
%  name        Which argument are we looking for?
%  value_type  What kind of value to return? (Default: 'bool')
%              'bool': True is returned if the name was found.
%              'val' : The next entry in the cell array is returned.
%                      Throws error if name was found last in cell array.
%              'num' : The next entry in the cell array is returned if it is numeric.
%                      Throws error if name was found last in cell array or followed by string.
%              'opt' : If the next entry in the cell array is not a string it is returned,
%                      otherwise default_if_no_following_value is returned.
%              default_if_not_named is always returned if the name is not found.
%              IMPROVEMENT: 'val*' and 'opt*' could return a cell array of
%              values, one per occurence of the name.
%  default_if_not_named          (Default: false) Used for 'val', 'num' & 'opt'.
%  default_if_no_following_value (Default: []) Used for 'opt'.
% RETURNS
%  value       See value_type for definition.
function value = get_argument_from_cells(options, name, value_type, default_if_not_named, default_if_no_following_value)
if nargin < 3
  value_type = 'bool';
end
if nargin < 4
  default_if_not_named = false;
end
if nargin < 5
  default_if_no_following_value = [];
end

matching_indices = strcmp(options, name);

switch value_type
  case 'bool'
    value = any(matching_indices);
    
  case {'val', 'num', 'opt'}
    index = find(matching_indices, 1);
    if any(index)
      value = default_if_no_following_value; % default if no value follows the name
      if length(options) > index
        if strcmp(value_type, 'val') % any value is used
          value = options{index+1}; % use the next array entry as value
        elseif strcmp(value_type, 'opt') % use if non-string, otherwise default
          if ~ischar(options{index+1}) 
            value = options{index+1}; % use the next array entry as value
          % else % use default
          end
        elseif strcmp(value_type, 'num') % require a number
          if isnumeric(options{index+1}) % number is always OK
            value = options{index+1}; % use the next array entry as number
          else
            error('get_argument_from_cells: No non-string value given for "%s" where value is mandatory.', name);
          end
        end
      elseif strcmp(value_type, 'val')
        error('get_argument_from_cells: No value given for "%s" (due to end of arguments) where value is mandatory.', name);
      end
    else
      value = default_if_not_named; % name not found
    end
end