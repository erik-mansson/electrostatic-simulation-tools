% chemical_fragment_possibilities(masses_to_identify, parent_formula)
%
% Wrapper to make chemical_identification() applicable to larger parent
% molecules than the current brute-force algorithm handles well,
% and to print result for multiple fragments by one call.
%
% The same list of fragment will be found, usually faster.
%
% In the printout ' e' is shown after fragment ions with an even number of electrons.
% Such cations (because the neutral fragment is a radical then) are favoured
% according to a heuristic rule found in chemistry texts/websites/Wikipedia. [Ref??]
%
% PARAMETERS
%  masses_to_identify   Mass number, or vector of mass numbers to identify
%  parent_masses        Chemical formula (string) of the parent molecule,
%                       or a row-vector with the mass number of each atom in it.
% EXAMPLES
%  chemical_fragment_possibilities(25, 'H_2 Na_2 O_2')
%  chemical_fragment_possibilities([15 19 38], 'C5N5H5')
%  chemical_fragment_possibilities(17, [1 1 16])
%
% SEE ALSO
%  chemical_identification The underlying method, used for one mass at a time
%  and without auto-reduction of parent molecule.
% Created by Erik P. Månsson 2015-09
%
function chemical_fragment_possibilities(masses_to_identify, parent_masses)

if ischar(parent_masses)
  parent_str = parent_masses;
  parent_masses = chemical_formula_to_mass(parent_masses);
else
  parent_str = strrep(chemical_formula_by_mass(parent_masses),'_',''); % skip underscores since this will be shown without TeX
end


% Count the number of each element in parent
[counts, elements] = hist(parent_masses,1:max(parent_masses));
mass_per_element = counts .* elements;

for i = 1:length(masses_to_identify)
  mass_to_identify = masses_to_identify(i);
  minimal_parent = parent_masses;
  
  % Reduce the count for this element so higher fragment than needed can't be produced
  % (using full parent may cause the brute-force search in chemical_identification() to be too slow).
  for e = find(mass_per_element > mass_to_identify)
    where = find(minimal_parent == e);
    minimal_parent(where(floor(mass_to_identify/e + 1):end)) = []; % delete the additional instances of element e
  end
  
  % TODO: exclude the hydrogens from parent and fragment, use a range of acceptable sub-masses that would
  parent_hydrogens = sum(minimal_parent == 1); % the number of hydrogens in parent
  if parent_hydrogens > 2 % (a threshold to not use this more complicated version when few hydrogens)
    dehydrogenated_parent = minimal_parent(minimal_parent > 1); % remove the hydrogens
    % For a mass in this range, adding hydrogens can reach mass_to_identify
    mass_range = [max(0,mass_to_identify-parent_hydrogens), mass_to_identify]; 
    %disp(sprintf('%d u (>= %d u without hydrogens) from de-hydrogenated minimal parent %s within %s:', mass_to_identify, mass_range(1), strrep(chemical_formula_by_mass(dehydrogenated_parent),'_',''), parent_str));
    disp(sprintf('%d u from de-hydrogenated minimal parent %s within %s:', mass_to_identify, strrep(chemical_formula_by_mass(dehydrogenated_parent),'_',''), parent_str));

%     tic
    chemical_identification(mass_range, dehydrogenated_parent);
%     disp(sprintf('%% Elapsed time (%ss) %ss per (dehydrogenated) binominal combination\n', texformat_SI(toc / 2.^length(dehydrogenated_parent), 2, false), texformat_SI(toc / 2.^length(minimal_parent), 2, false) ));

  else % Version 1:
    
    if length(minimal_parent) == length(parent_masses)
      disp(sprintf('%d u from parent %s:', mass_to_identify, parent_str));
    else
      disp(sprintf('%d u from minimal parent %s within %s:', mass_to_identify, strrep(chemical_formula_by_mass(minimal_parent),'_',''), parent_str));
    end

    % tic
    chemical_identification(mass_to_identify, minimal_parent);
    % disp(sprintf('%% Elapsed time %ss per binominal combination\n', texformat_SI(toc / 2.^length(minimal_parent), 2, false)) );

  end
  
end
