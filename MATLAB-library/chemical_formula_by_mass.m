% formula = chemical_formula_by_mass(masses, group_by_element)
%
% Returns a readable form for the chemical formula of a molecule,
% from input giving the MASS number (not nuclear charge) of its atoms.
%
% IMPROVEMENT: is there some rules for the order of elements,
% e.g. why is the conventional names CO_2 and H_2O not of the same form??
%
% PARAMETERS
%  masses
%  group_by_element (Default: true)
%                   when true:  [16 12 16] becomes O_2C
%                   when false: [16 12 16] becomes OCO
% SEE ALSO
%   chemical_formula_to_mass
function formula = chemical_formula_by_mass(masses, group_by_element)

% Known elements by MASS number (Hydrogen isotopes included)
elements = {'H','D','T','He','?',...
    'Li','^7Li','?','Be','^{10}B','B','C','^{13}C',...
    'N','^{15}N','O','^{17}O','^{18}O','F','Ne','^{21}Ne','^{22}Ne'};
elements{23}='Na';
elements{24}='Mg';
elements{27}='Al';
elements{28}='Si';
elements{32}='S';
elements{35}='^{35}Cl';
elements{36}='^{36}Ar';
elements{37}='^{37}Cl';
elements{38}='^{38}Ar';
elements{40}='Ar';

if nargin == 0 && nargout == 1
  % Special case: if no arguments are given, return the list of known elements!
  % For use by chemical_formula_to_mass.m
  formula = elements;
  return;
end
if ischar(masses) % if a string was given, just return it
  formula = masses;
  return;
end

masses(masses == 0) = []; % remove zeroes, if they would occur

if nargin < 2
  group_by_element = true;
end
formula = '';
used_masses = false(1,max(masses));
for m = masses
  if ~used_masses(m)
    if group_by_element
      count = sum(masses == m);
      used_masses(m) = true;
    else
      count = 1;
    end
    if m <= length(elements)
      elem = elements{m};
    else
      error('The chemical element with M=%d is not named in the elements array', m);
    end
    if count > 1
      formula = [formula elem '_' num2str(sum(masses == m))];
    else
      formula = [formula elem];
    end
  end
end
