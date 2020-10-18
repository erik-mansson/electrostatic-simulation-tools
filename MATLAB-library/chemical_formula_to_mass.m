% [masses, molecular_weight] = chemical_formula_to_mass(formula)
% Convert a chemical formula to an array of integer mass numbers, one element per atom.
%
% PARAMETERS
%  formula    String with chemical symbols and numeric suffixes with or without underscores.
%             Note that only the first character in an element is uppercase.
%             Isotope prefixes need to include the curly braces, e.g. '^{13}C'.
%             Parentheses are currently not supported.
%             Examples: 'OH', 'H2O', 'H_2O', 'H2Na2O2', 'H_2 Na_2 O_2', 'Al2^{13}CB'.
% RETURNS
%  masses            Row vector with mass number of each chemical element,
%                    following the order and grouping in the formula.
%  molecular_weight  The same as sum(mases)
%
% SEE ALSO
%   chemical_formula_by_mass
function [masses, molecular_weight] = chemical_formula_to_mass(formula)

% Remove meaningless but allowed charactes
%s = regexprep(formula, '[0+- _]', '');

% Get list of groups, where each group is the longest possibility
% of (Element)_(number) or (Element)(number) or (Element)
% Note that only the first character in an element is uppercase.
% [groups, remaining_text] = regexp(formula, '([A-Z][a-z]*)(_?\d+)?', 'tokens', 'split'); % without isotope prefixes
[groups, remaining_text] = regexp(formula, '(([^][{]\d+[}])?[A-Z][a-z]*)(_?\d+)?', 'tokens', 'split');

for i = 1:length(remaining_text)
  if ~isempty(strtrim(remaining_text{i}))
    error('The part "%s" in "%s" does not look like a chemical element (e.g. Na B O Uut).', remaining_text{i}, formula);
  end
end

% Get cell array of elements/isotopes
element_list = chemical_formula_by_mass();
% Cache searches so the same string won't be searched for again
element_mass = containers.Map('KeyType','char', 'ValueType','int32');

masses = [];
for i = 1:length(groups)
  element = groups{i}{1};
  
  % Get the mass
  if element_mass.isKey(element)
    mass = element_mass(element);
  else
    % Search the list
    mass = find(strcmp(element_list, element)); % index is equal to mass number, currently not supporting multiple elements at same mass!
    if isempty(mass)
      error('The chemical element "%s" in the formula "%s" is not supported. Correct the formula or add this element in chemical_formula_to_mass.m.', element, formula);
    else % OK
      % Save to cache to not search array again for this element
      element_mass(element) = mass;
    end
  end
  
  % Convert the count
  if isempty(groups{i}{2})
    count = 1;
  else
    if groups{i}{2}(1) == '_'
      % Skip the underscore
      groups{i}{2} = groups{i}{2}(2:end);
    end
    count = str2double(groups{i}{2});
    if isnan(count) || count < 1 || (round(count) ~= count)
      error('The part "%s" (of "%s%s" in "%s") appears to be an invalid number. Only positive integers can be used.', groups{i}{2}, element, groups{i}{2}, formula);
    end
  end
  
  % Append to output
  masses = [masses, repmat(mass,1,count)];
end

molecular_weight = sum(masses);
