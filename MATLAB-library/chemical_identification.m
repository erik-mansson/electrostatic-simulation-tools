% List molecules with a given molecular weight (integer mass number)
% that can be formed as a fragment of the given parent molecule.
%
% In the printout ' e' is shown after fragment ions with an even number of electrons.
% Such cations (because the neutral fragment is a radical then) are favoured
% according to a heuristic rule found in chemistry texts/websites/Wikipedia. [Ref??]
%
% PARAMETERS
%  mass              integer mass to identify
%                    Optionally as an optimization for hydrogen-containing molecules
%                    a two-element array [m-h, m] can be given. Then any fragment
%                    whose mass f satisfies m-h <= f <=m will be taken as a mach.
%                    At most h hydrogen atoms will be appended to it to reach the mass m.
%
%                    given in case 
%  parent            string  chemical formula or array of mass numbers for the atoms.
%  group_by_element  boolean (default: true)
%                    Use numbers in output chemical formulas, e.g. C_2H_2 rather than CCHH.
%
% RETURNS
%  candidate_formulas   N-by-1 cell array of chemical formulas for the found molecular candidates.
%                       If underscores are not desired (i.e. not printing/drawing using TeX formatting)
%                       they can be removed by the expression strrep(candidate_formulas, '_','').
%  candiate_masses      N-by-1 cell array of row-vectors, where each row-vector
%                       contains the mass numbers of elements included.
%
% SEE ALSO
% chemical_formula_to_mass chemical_formula_by_mass
% chemical_fragment_possibilities
%                    For wrapper to call this function in a more efficient way for moderately large molecules.
%                    The same list of fragment will be found, often faster.
% Created by Erik P. Månsson 2015-09
%
function [candidate_formulas, candiate_masses] = chemical_identification(mass, parent, group_by_element)
if nargin < 3
  group_by_element = true;
end

if size(parent,1) ~= 1 || isempty(parent)
  error('The parent mass array should be a row-array, or a string with the chemical formula.');
end
if ischar(parent)
  parent_masses = chemical_formula_to_mass(parent);
elseif isnumeric(parent) && all(parent == round(parent))
  parent_masses = round(parent);
end
if isempty(mass) || size(mass,1) ~= 1 || size(mass,2) > 2
  error('The mass to search for should be a scalar, or an 1-by-2 array of [min mass without hydrogens, desired mass].');
end
if any(mass ~= round(mass))
  error('The molecular weight (mass number) must be an integer.');
end
max_mass = mass(end);
if length(mass) == 2
  min_mass = mass(1); % min mass without hydrogens
  if any(parent_masses == 1)
    error(sprintf('When optimized hydrogen treatment is used, i.e. searching for molecule without hydrogens then adding \nas many as needed, the parent molecule formula should not include any hydrogens.'));
  end
else % searching for a single mass, no special treatment of hydrogen atoms
  min_mass = max_mass;
end
if max_mass < 1
 error('The molecular weight (mass) to find must be a positive integer.');
end
if min_mass < 0
  min_mass = 0; % Can't use all hydrogens, because target mass is lower than number of hydrogen atoms in parent
end
if min_mass > max_mass
  error('The min mass without hydrogens may not be larger than the desired mass.');
end
N = length(parent_masses);

% Brute-force approach. Try all combinations and append those that match the given mass.
% By treating the include-array as a binary number, and iterating in "numerical" order
% (now on the reversed parent_masses) all combinations can be tried.
% indices = N:-1:1; % to get all bits 
indices = 1:N; % to get all bits. This order makes iteration happen as if a "reversed" parent_masses was used, but it is OK and it seems plausible that bitget() may have been optimized for this order.

if N <= 32
  include_bitmask = uint32(0);
  end_value = uint32(2^N - 1);
  parent_masses= uint32(parent_masses); % for multiplication with the include-variable to work
elseif N <= 64
  include_bitmask = uint64(0);
  end_value = uint64(2^N - 1);
  parent_masses= uint64(parent_masses); % for multiplication with the include-variable to work
else
  error('The current implementation can not handle more than 64 elements in the parent molecule.\n As a workaround You could use fake elements like 18 for H_2O if that part is expected to remain unfragmented.');
end

candidate_formulas = {};
candiate_masses = {};
print_results = nargout == 0;

while true %include_bitmask <= end_value
  
  include = bitget(include_bitmask, indices); % get the relveant bit
  
  weight = sum(parent_masses .* include);
  if weight >= min_mass && weight <= max_mass
    % Matching candidate found
    candidate = parent_masses(include == 1);
    
    % Append hydrogens if allowed and needed
%     extra_hydrogens = max_mass - weight;
%     if extra_hydrogens > 0
    candidate = [candidate, ones(1, max_mass - weight)];
    
    % Version 1: logging by formula. Although conversion to formula can be expected to be slow,
    % this part is not run very often (the number of isomers per new fragment is not disturbingly high).
    formula = chemical_formula_by_mass(candidate, group_by_element);
    if ~any(strcmp(candidate_formulas, formula))
      % The candidate has not been found before
      candidate_formulas{end+1,1} = formula;

%     % Version 2 (was slightly slower): logging by mass-array (which could print isomers of same fragment if parent molecule was given in a structured form (not just sum formula))
%     % Check whether the candidate is new (has not been found before)
%     % search cell array of mass vectors
%     new = true;
%     for check_index = 1:numel(candiate_masses)
%       if isequal(candiate_masses{check_index,1}, candidate)
%         new = false; break;
%       end
%     end
%     if new
%       candiate_masses{end+1,1} = candidate;
%       formula = chemical_formula_by_mass(candidate, group_by_element);
%       candidate_formulas{end+1,1} = formula;

      if print_results
        % Print results without underscores, to ease reading when TeX formatting is not available.
%         disp(strrep(formula, '_',''));
        
        candidate(candidate == 12 | candidate == 16) = 0; % carbon and oxygen has even number of electrons
        candidate(candidate == 14) = 1; % Nitrogen (and hydrogen) have odd number of electrons
        if all(candidate <= 1) % only HCNO atoms, then can print info from even/odd electron count rule
          electron_count = sum(candidate) - 1; % subtract one since cation considered
          if mod(electron_count, 2) == 0
            % Even number of electrons. Such cations are favoured according to the rule (radical neutral fragment instead).
            disp([strrep(formula, '_','') ' e']);
          else
            % Odd number of electrons.
%             disp([strrep(formula, '_','') ' -o']);
            disp(strrep(formula, '_',''));
          end
        else
          disp(strrep(formula, '_',''));
        end
      else
        if nargout >= 2
          candiate_masses{end+1,1} = candidate;
        end
      end
    end
  end
  
  if include_bitmask >= end_value
    % In case end_value is also the maximum value for the data type used (int32 or int64),
    % the incremented include_bitmask will by Matlab's overflow handling stay at its current value.
    % Therefore this extra check is needed to ensure the loop won't run forever.
    break;
  end
  % Increment the inclusion mask
  include_bitmask = include_bitmask + 1;
end

if nargout == 0
  
  % Called (from command line) with semicolon after, not getting output variable.
  if isempty(candidate_formulas)
    disp(sprintf('No fragment with mass number %d wcould be produced from the parent molecule.', max_mass));
  else
    % Print results without underscores, to ease reading when TeX formatting is not available.
    % disp(strrep(candidate_formulas, '_',''));
%     shortened = strrep(candidate_formulas, '_','');
%     for i = 1:length(candidate_formulas)
%       disp(shortened{i,1}); % to avoid quotes around each formula
%     end
  end
end
