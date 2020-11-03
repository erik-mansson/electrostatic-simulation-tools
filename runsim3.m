% [flys, logfile_path, DLT_path] = runsim3(workbench_name, mass_in_u, energies, source_point_pattern, variables, runsim3_job_index)
%
% New approach to running SIMION:
%  * A .fly2 file is created based on a template and parameters.
%    Thanks to improved workbench Lua script it can loop over groups
%    of particles with different properties (currently only implementing for kinetic energies).
%  * An adjust.lua can be created for any other overrides or Lua-code given after variables{'script'}
%    (This will need to be the same for all parallell instances.)
%  * SIMION is run without the commandline2.lua
%
% A global variable runsim3_counter is also used to ensure unique temporary filenames.
% It may be reset to [] before starting a batch job, if it is certain that no other job is running in parallel.
%
% PARAMETERS
%   workbench_name        Name of SIMION workbench, with or without the .iob suffix.
%                         Currently assumed to be within predetermined simulation_root directory.
%   mass_in_u             [u] mass of particle. Charge is by default sign(mass_in_u), but if |mass| < 1u is specified,
%                         then the electron mass will be used and the charge set to negative.
%   energies              [eV] array of kinetic energies to simulate for.
%   source_point_pattern  If a single integer index >= 0 is given: Use a predefined source point pattern
%                         (see the workbench Lua script for available pattern).
%                         If an array [std_axial std_transverse] or [std_axial std_transverse std_transverse_other]
%                         is specified: use a pseudorandom (reproducible) 3D gaussian with these standard deviations
%                         along the SIMION X,Y,Z coordinates. If std_transverse_other is omitted it defaults tostd_transverse.
%                         If an empty array is given, the default is used.
%                         Default: currently [0.06 1 0.06] 3D gaussian.
%   variables             Cell array of name,value pairs with additional parameters.
%
%                         Simulation settings (default):
%                         '_number_of_directions' (8)
%                         'particles_per_direction' (4)
%                         'debug' (false) -- retain the job files afterwards instead of deleting them
%                         'trajectories' (dependent on workbench .rec, often false) 
%                             tells whether to record individual particle start & impact coordinates.
%                         'initial_vy' (0) constant offset in intial velocities [mm/us]=[km/s] e.g. to represent speed of molecular beam
%                         'initial_v_vertical' (0) constant offset in intial velocities [mm/us]=[km/s]
%                             e.g. to represent speed of molecular beam.
%                             NOTE: For consistencity with coordinate swapping 'initial_v_vertical' 
%                             is recommended here as an alias for 'initial_vz' which actually
%                             refers to vertical (up, usually polarization axis)
%                             while z in most Matlab variables (e.g. source_point_z) 
%                             refers to horizontal (TOF axis).
%
%                         '_isotropic' (default depends on workbench)
%                               0: discrete directions 0 to 360 degrees
%                               1: isotropic, random distribution
%                               2: 2D-VMI-style directions +90 to -90 degrees 
%                                  (combining particles with z-momentum in forward and backward directions)
%                               3: 3D-VMI-style directions +90 to -90 degrees 
%                                  (combining particles with z-momentum in forward and backward directions),
%                                  (discretized,approx.) isotropic solid angle density by varying number of particles per ray
%                               4: 3D-VMI-style directions not including p_Z=0 (+-90 degrees)
%                                  (combining particles with z-momentum in forward and backward directions),
%                                  (discretized,approx.) isotropic solid angle density by varying number of particles per ray
%
%                         Potentials via parameters: '_U_Bfirst', '_U_Afirst', ...
%                         'potentials': Fully customized potentials, given as a row-array following the SIMION indexing of electrodes.
%
%                         'script' Other custom Lua script whose source code will be written to adjust.lua and performed
%                         during initialization of the workbench script. NOTE: the 'script' variable is not thread-safe, 
%                         meaning you must avoid running parallel calls to runsim3 if it is used.
%
%   job_index             (Default: 0) If parallell loops are used, this number should be different from all
%                         concurrently running instances to ensure no filename collisions occur.
%
% RETURNS
%  flys                   Structure with the information parsed from the logfile.
%  logfile_path           Absolute path to the logfile.
%  DLT_path               (TODO) Absolute path to the DLT file. Requesting this output causes DLT writing to be enabled.
%
% EXAMPLE
%  runsim3('example', 0, [1 3])
%
% SEE ALSO
%  readsim3.m, script_VMI_sim_and_calib.m, optimize_potentials_objective.m
function [flys, logfile_path, DLT_path] = runsim3(workbench_name, mass_in_u, energies, source_point_pattern, variables, runsim3_job_index)

% Strip trailing .iob (if given) from workbench name
workbench_name = regexprep(workbench_name, '\.iob$', '');

simulation_root = '.\example model\';  % Assume that the model files are in a path relative to current directory (not usually the case, but better than defaulting to an invalid path)
% simulation_root = 'C:\Users\<USERNAME>\Documents\SIMION-models\'; % TODO: To be customized to relevant path for each user
simion_path = 'C:\Program Files\SIMION-8.1\simion.exe'; % TODO: To be customized to relevant path for each user


% Simulate only singly charged particles, using +1 or -1 depending on sign & magnitude of mass.
charge_in_e = sign(mass_in_u);
mass_in_u   = abs(mass_in_u); % restore to positive mass for all
which_electrons = abs(mass_in_u) < 1; % Zero or small mass is used as a shortcut for electron mass
mass_in_u(which_electrons) = 0.000548579903; % Using SIMION's default electron mass, wich to eight digits agrees with later CODATA value.
charge_in_e(which_electrons) = -1;
% IMPROVEMENT: maybe allow other charge via variables{} or second row of mass_in_u?

if nargin < 3
  energies = [0.4, 1.95, 3.5, 5.05, 6.6]; % [eV] 1.55 eV spacing
end
energies = energies(:)'; % ensure shaped as row vector
if length(energies) ~= 1 && length(mass_in_u) ~= 1 && length(energies) ~= size(mass_in_u,2)
  error('If multiple energies and multiple masses are given, the array sizes must match. Got %d energies and %d masses.', length(energies), size(mass_in_u,2));
end
  
if nargin < 4 || isempty(source_point_pattern) || all(source_point_pattern < 0)
  source_point_pattern = [0.06 1 0.06]; %[mm] Default: 3D-Gaussian via fly2
end
if nargin < 5
  variables = {}; 
end
if nargin < 6
  runsim3_job_index = 0;
end

% Provide default source_point_z for the listed workbenches:
source_point_y = 0; %[mm] default to which source_point_offset_vertical is added. DEBUG: before 2020-02-25 the offset was absolute from 0, not relative to model's default.
% Previously source_point_offset_vertical had to be set to a geometry-dependent value, the geometry-offset was not added here. But makes more sense to let "0" mean "nominal central point for geometry"
source_point_z = [];
source_point_vertical = 0; %[mm] default to which source_point_offset_vertical is added
required_number_of_potentials = 0; % default no info about how many potentials, accept anything

switch workbench_name
  % It is recommended (required for some usage) to give some geometrical or electrode-array 
  % info about each model here even if it repeats info available in the Lua script of the model.
  case {'example'}
    source_point_z = 8.5; % [mm] midpoint between repeller and extractor in simulation geometry
    required_number_of_potentials = 3;
  case {'ecyl8e'}
    source_point_z = 640-550; % [mm]
    required_number_of_potentials = 10;
  case {'both8eu', 'both8eu_B2'}
    if strcmp(workbench_name, 'both8eu_B2')
      required_number_of_potentials = 19;
    else
      required_number_of_potentials = 18;
    end
    source_point_z = 640; % [mm]
    source_point_y = 65; % [mm]
    source_point_vertical = 70; % [mm] default for this geometry (which doesn't vertical mirror-symmetry)
  
  otherwise
end
normal_source_point_z = source_point_z;

% Override source_point_z with variable, if given. An empty array [] remains if also no default was given for the workebench_name.
% For 3D-gaussian source_point_pattern, [] will raise an error, while for other patterns the absence of setting the source_point_z adjustable
% means that the default defined in the Lua script will be used.
source_point_z = get_argument_from_cells(variables, 'source_point_z', 'num', normal_source_point_z, normal_source_point_z); 
if get_argument_from_cells(variables, 'source_point_offset_z', 'num')
  if isempty(source_point_z)
    error('source_point_offset_z can not be used if source_point_z is not defined for the model name "%s".', workbench_name);
  end
  if get_argument_from_cells(variables, 'source_point_z', 'num')
    warning('source_point_offset_z and source_point_z should not be given simultaneously, let source_point_z be given by workbench_name.');
  end
  % If source_point_z is known based on model, and 'source_point_offset_z' is given,
  % and no 'source_point_z' override is given among variables,
  % then we use it as a user-friendly way of modifying the source_point_z
  % (Lua or SIMION never handle the relative source_point_offset_z, only the absolute source_point_z).
  offset = get_argument_from_cells(variables, 'source_point_offset_z', 'num', [], []);
  if isfinite(offset)
    fprintf('Using source_point_offset_z to shift source_point_z from %.0f to %.0f.\n', source_point_z, source_point_z + offset);
    source_point_z = source_point_z + offset;
    % Remove the offset variable, since we no longer want to send offsets to Lua
    index = find(strcmp(variables,'source_point_offset_z'));
    variables(index:index+1) = []; % remove used variable that should not be sent to Lua
  end
end


absolute = get_argument_from_cells(variables, 'source_point_y', 'num', [], []);
offset = get_argument_from_cells(variables, 'source_point_offset_y', 'num', [], []);
need_to_specify_y = false;
if ~isempty(absolute)
  % In case an absolute coordinate was given:
  if ~isempty(offset)
    error('Both an absolute source_point_y=%.1f mm and relative source_point_offset_y=%.1f mm were given. From 2020 an offset is converted to absolute by runsim3 by adding the default for each geometry (=%.1f mm for %s).', absolute, offset, source_point_y, workbench_name);
  else
    source_point_y = absolute;
    % no need to append it to variables, it is already there (thus not setting need_to_specify_y to true)
  end
else
  if ~isempty(offset)
    if offset ~= 0
      % In case relative offset was given (NOTE: before 2020-02-25 the _offset variables were actually absolute, but it seems more user-friendly to let offset=0 mean nominal point for any geometry)
      if abs(offset) >= 15 || (~isempty(absolute) && absolute == offset)
        warning('runsim:y_offset_relative', 'Suspicously large source_point_y_offset=%.1f mm in model %s with default source_point_y=%.1f mm. 2020-02-25 runsim3 was changed so that y_offset is relative to the geometry''s default, not replacing it.', offset, workbench_name, source_point_y)
      end
      source_point_y = source_point_y + offset; % realies on the default source_point_y having been defined for the geometry in the switch-case above
      need_to_specify_y = true;
    end
    % Remove the offset variable, since we no longer want to send offsets to Lua (and the variable is renamed to source_point_y in newer scripts)
    index = find(strcmp(variables,'source_point_offset_y'));
    variables(index:index+1) = []; % remove used variable that should not be sent to Lua
  end
end
  
% Since the x and z are confusingly swapped sometimes and sometimes not, 
% I now use the name "vertical" for the SIMION Z (Matlab x) which is 
% the coordinate perpendicular to both optical axis and TOF-axis.
absolute  = get_argument_from_cells(variables, 'source_point_vertical', 'num', [], []);
offset = get_argument_from_cells(variables, 'source_point_offset_vertical', 'num', [], []);
need_to_specify_vertical = false;
if ~isempty(absolute)
  % In case an absolute coordinate was given:
  if ~isempty(offset)
    error('Both an absolute source_point_vertical=%.1f mm and relative source_point_offset_vertical=%.1f mm were given. An offset is converted to absolute by runsim3 by adding the default for each geometry (=%.1f mm for %s).', absolute, offset, source_point_vertical, workbench_name);
  else
    source_point_vertical = absolute;
    % no need to append it to variables, it is already there (thus not setting need_to_specify_vertical to true)
  end
  index = find(strcmp(variables,'source_point_vertical'));
  variables{index:index+1} = []; % remove used variable that should not be sent to Lua
elseif ~isempty(offset)
  if offset ~= 0
    source_point_vertical = source_point_vertical + offset; % realies on the default source_point_y having been defined for the geometry in the switch-case above
    need_to_specify_vertical = true;
  end
  % Currently the Lua scripts do not handle vertical offset, i.e. it only works with 3D-Gaussian source point pattern for which a FLY2 file is created.
  % -- otherwise if need_to_specify_vertical: variables{end+1} = 'source_point_vertical'; variables{end+1} = source_point_vertical;
  index = find(strcmp(variables,'source_point_offset_vertical'));
  variables(index:index+1) = []; % remove used variable that should not be sent to Lua
end


if length(source_point_pattern) >= 2
  % 3D-gaussian volume, using SIMION's generator rather than generator in Lua script.
  % SIMION's pseudorandom generator is independent of the adjustable random seed variable
  % and will not be reproduceable between groups (energies) within one simulation, but the simulation as a whole is reproducable.
  % I.e. if a run contains two identical groups their random start position will differ,
  % but running the two-group simulation again doesn't change it further.
  std_axial = source_point_pattern(1); %[mm] distances should be given in millimetres
  std_y = source_point_pattern(2); %[mm] distances should be given in millimetres, for vertical axis in simulation
  if length(source_point_pattern) == 2
    std_x = std_y;
  elseif length(source_point_pattern) >= 3
    std_x = source_point_pattern(3); %[mm] distances should be given in millimetres
  else
    error('At most three dimensions of source_point_pattern is supported.');
  end
  
  index = find(strcmp(variables,'source_point_z')) + 1;
  if ~isempty(index)
    fprintf('Using source_point_z variable of %.6g mm instead of %.6g.\n', variables{index}, source_point_z);
    source_point_z = variables{index};
    % Since source_point_z is set in the .fly2 file, no need to pass it on to the Lua script (default there will have no effect).
    variables(index-1:index) = [];
  elseif isempty(source_point_z) % had no default for the geometry (workebench_name)
    error('The source_point_z variable [mm] must be specified when using Gaussian source volume. You may add a default for the workbench "%s" by editing runsim3.m or add a ''source_point_z'',value pair in variables.', workbench_name);
  end
  
  need_to_specify_y = true;

  position_str = sprintf('gaussian3d_distribution {mean = vector(%.3f, %.3f, %.3f), stdev = vector(%.4f, %.4f, %.4f)}', ...
                          source_point_z, source_point_y, source_point_vertical, ...
                          std_axial, std_y, std_x);
  position_description = sprintf('gaussian std [%.5g,%.5g,%.5g]', std_axial, std_y, std_x);
  source_point_pattern = -1; % to disable position setting by workbench Lua script
  
  % Indicate that _source_point_spacing wasn't used
  index = find(strcmp(variables,'_source_point_spacing')) + 1;
  if ~isempty(index)
    variables{index} = -1;
  else
    variables{end+1} = '_source_point_spacing'; variables{end+1} = -1;
  end

else
  % Deterministic source points or pseudorandom by generator in Lua script.
  % If pseudorandom, the seed is set to make each group reproducable and independent of others.
  position_str = 'vector(0, 0, 0) --  Will be overridden by Lua script anyway'; % since source_point_pattern is not a 3D-gaussian
  position_description = '';
  
  % In this case the Lua adjustables for position are not always assigned from here, to allow the script's default to be used.
  % 'source_point_z' is listed in allowed_adjustables and will be passed on without change, if given.
  % For the vertical and horizontal we allow the user to give only an offset with respect to the nominal point for the geometry,
  % see above.
  
  if need_to_specify_vertical
    error('Currently the Lua scripts do not handle vertical offset, i.e. it only works with 3D-Gaussian source point pattern for which a FLY2 file is created.')
  end
end

if need_to_specify_y
  % Especially for non-3D-Gaussian source_point_pattern, we need to pass it as adjustable variable to Lua, but OK to always do it.
  % Also needs to be done in case changed by a nonzero source_point_offset_y
  variables{end+1} = 'source_point_y'; variables{end+1} = source_point_y;
  % NOTE: this will raise error in old Lua scripts where the adjustable was called "source_point_offset_y" even if it was absolute. Just rename the adjustable in such scripts!
end

% (Now happens by _source_point_pattern not being in list of allowed_adjustables:)
% % Ensure that the _source_point_pattern in Lua is controlled by the source_point_pattern argument here.
% index = find(strcmp(variables,'_source_point_pattern')) + 1;
% if ~isempty(index)
%   variables{index} = source_point_pattern;
% else
%   variables{end+1} = '_source_point_pattern'; variables{end+1} = source_point_pattern;
% end

index = find(strcmp(variables,'_number_of_directions')) + 1;
if ~isempty(index)
  number_of_directions = variables{index};
else % defaults
  number_of_directions = 8;
end
index = find(strcmp(variables,'particles_per_direction')) + 1;
if ~isempty(index)
  particles_per_direction = variables{index};
  variables(index-1:index) = []; % this variable is not to be sent to the Lua script
else
%   particles_per_direction = 50;
  particles_per_direction = 4; % DEBUG something small while testing
end

fly2 = sprintf(['_G.FLY2_IONS_PER_GROUP = %d\n_G.FLY2_SOURCE_DISTR = "%s"\n' ... % Let the workbench Lua script know after how many ions a new energy starts. And know gaussian standard deviations.
'local t = {coordinates = 0}\n'], number_of_directions*particles_per_direction, position_description);
used = false(size(mass_in_u));
%'local energies = {%s}\n' ...

for mi = 1:length(mass_in_u)
  if length(mass_in_u) == 1 || length(energies) == 1
    % Only one mass or one energy (repeat it for all energies|masses)
    which = true(size(energies));
  elseif length(mass_in_u) ~= length(energies)
    error('If multiple masses are and multiple energies are given, it must be as two arrays of the same size. If a single mass is given, it will be used for all energies.\n masses: %s, energies: %s', mat2str(mass_in_u), mat2str(energies));
  else
    % Arrays of same length for mass and energy. Make one particle group per mass,energy combination.
    % In Lua FLY2 this is further grouped to one command block per mass (listing all its energies):
    which = ~used & mass_in_u==mass_in_u(mi) & charge_in_e==charge_in_e(mi);
    if ~any(which)
      % just another energy for the same mass already used
      continue;
    end
  end
  % Use a loop to build the table passed to "particles", with one entry per kinetic energy. One loop per (mass,charge)-combination.
  fly2 = sprintf(['%sfor index,energy in ipairs({%s}) do\n' ... 
    '  t[#t+1] = standard_beam {n = _G.FLY2_IONS_PER_GROUP, tob = 0, cwf = 1,\n' ...
    '    mass = %.12g, charge = %d, ke = energy,\n' ...
    '    color = (index <= 2) and index-1 or index,\n' ... % Means: index 1 ==> color 0, 2==>1, otherwise color = index (not using color 2 which is for equipotential curves)
    '    direction = vector(1,0,0), -- Will be overridden by Lua program\n' ...
    '    position = %s\n  }\nend\n'], ...
    fly2, sprintf('%.4g,', energies(which)), mass_in_u(mi), charge_in_e(mi), position_str);
  used(which) = true; % don't reuse entries
end
fly2 = sprintf('%sparticles(t)\n', fly2);

debug_mode = false; % can be set to true by 'debug' parameter

% (Skipped, globals not good within parallel for loop)
% % Define a unique (within Matlab session) number for the SIMION launch, to possibly run several in parallell
% global runsim3_job_index;
% if isempty(runsim3_job_index)
%   runsim3_job_index = 0;
% end

job_name = sprintf('sim_%05d%s', runsim3_job_index, char(floor(65 + 25*rand(1,3))) );
fly2_name = [job_name '.fly2'];
output_name = [job_name '.txt'];
% (when global:) runsim3_job_index = runsim3_job_index + 1;

further_arguments = '';
% TODO: prepare adjustment.lua to set parameters for potentials (do its writing near executing the command)
% IMPROVEMENT: If desirable to run instances in parallell, then adjustment.lua name needs to be variable. Either via another global assigned from the .fly2 or by using a variably-named workbench.lua, setting the variables of interest and then perhaps able to load the main workbench Lua file (or just write altered version of big script for every run?)
% IMPROVEMENT: Figure out if it would be possible & better to use the SIMION fastadj command instead of the adjustable parameters. It would mean Matlab scripts need to know more details of spectrometers, but would allow more fine-grained control (though probably not needed for a VMIS)

% NOTE/IMPROVEMENT: Can set completely custom potentials like this in adjustment.lua (would something in .fly2 work too, to reduce number of files written? Maybe using a global U to message?):
%  adjustable DoNotSetPotentials = -1
%  U = {0,  0,nil,nil,nil, 0, 0,-90,-73.33,-56.67,-40, -62,-62,-62,0,12,24,36,48,60,72,84,96,108,120,120}
% But currently concentrating on the few-parameter adjustment via Lua workbench program.

allowed_adjustables = { ...
  '_U_Bfirst','_U_Afirst', '_U_Bfree', '_U_Afree', '_U_Blast', '_U_Alast', '_A_free', '_A_plateau', '_U_BMCP', '_U_AMCP', ...
  '_number_of_directions', 'source_point_z', 'source_point_y', '_source_point_spacing', '_isotropic', ...
  'initial_vy', 'initial_v_vertical', 'multipole_order' % New 2019-06, so far only in ecyl8e but can be pasted in other models too
  };
% NOTE: only late geometries support initial_v_vertical -- in earlier geometries it would be ignored or give error (didn't test)
% NOTE: the naming was a bit confused. 'initial_vz' refers to vertical so it has been renamed to 'initial_v_vertical' 
% (ought to be renmaed to initial_v_vertical) but source_point_z refers to horizontal.
if any(strcmp(variables, 'initial_vz'))
  % Handle old code by just renaming the variable.
  variables{strcmp(variables, 'initial_vz')} = 'initial_v_vertical';
end
initial_v_vertical = get_argument_from_cells(variables, 'initial_v_vertical', 'num', 0); % [m/s] store this value which we wil put in the output source.initial_v_vertical

i = 1;
adjustment_script = '';
while i < length(variables) % (stop before i=length(variables) as index after is needed for value)
  name = variables{i};
  OK = any(strcmp(allowed_adjustables, name));
  if ~OK && name(1) ~= '_'
    OK = any(strcmp(allowed_adjustables, ['_' name]));
    if OK % Prepend _ to name (used in SIMION for adjustables that are shown in GUI)
      name = ['_' name];
    end
  end
  if OK
    % Append adjustment to command line
    further_arguments = sprintf('%s --adjustable %s=%.8g', further_arguments, name, variables{i + 1});
    % Remove name,value pair from list of unprocessed variables
    variables(i:i+1) = [];
  else
    switch name
      case 'script'
        adjustment_script = [adjustment_script variables{i + 1}];
        if ~isempty(adjustment_script) && runsim3_job_index ~= 0
          error('The script variable for custom adjustment is currently not thread-safe. Calling runsim3 with a job_index %d suggests that multiple instances are being launched in parallel.', runsim3_job_index)
        end
        variables(i:i+1) = []; % Remove name,value pair from list of unprocessed variables
        
      case 'potentials'
        % Set the array of potentials directly, not via the (smaller) set of _U_...-parameters.
        % This allows more generic SIMION geometries to be used, without further complexity in the Lua script.
        if isempty(variables{i + 1}) || ischar(variables{i + 1}) || any(~isfinite(variables{i + 1}))
          error('Got invalid potentials: %s for %s (%s).', num2str(variables{i + 1}(:)', '%.8g, '), workbench_name, job_name);
        end
        if required_number_of_potentials ~= 0 && length(variables{i + 1}(:)) ~= required_number_of_potentials
          error('Got %d potentials but %d are required for the workbench %s', length(variables{i + 1}(:)), required_number_of_potentials, workbench_name);
        end
        % Using the .fly2-script here rather than adjust.lua, to be thread-safe (unique name).
        % This works via global _G.U instead of just U, with support from the workbench Lua script.
        % The U[0] entry is just needed when workbench Lua script has some electrode index (e.g. B_MCP) set to 0 for (meaning there is no B MCP electrode in this geometry)
        fly2 = sprintf('%s_G.U = {%s}; _G.U[0]=0;\n', fly2, num2str(variables{i + 1}(:)', '%.8g, '));
        % Prevent other adjustables from overrideing the potential array
        further_arguments = sprintf('%s --adjustable DoNotSetPotentials=-1', further_arguments);

        variables(i:i+1) = []; % Remove name,value pair from list of unprocessed variables
        
      case 'debug'
        debug_mode = any(variables{i + 1});
        % In debug mode, temporary files are not deleted
        variables(i:i+1) = []; % Remove name,value pair from list of unprocessed variables
        
      case 'trajectories'
        if any(variables{i}) % if value is true
          % Enable recording of the individual particle trajectories
          %further_arguments = sprintf('--recording=record.rec %s', further_arguments); % old until 2015-12-29, did not record third coordinate (z in SIMION) which is zero unless _isotropic==1 or >=3
          further_arguments = sprintf('--recording=record3.rec %s', further_arguments); % old until 2015-12-29, did not record third coordinate (z in SIMION) which is zero unless _isotropic==1 or >=3
        else % 'trajectories',false given as key,value-pair
          % Ensure no trajectories get recorded. Disable non-Lua output by specifying nonexistent .REC-file
          further_arguments = sprintf('--recording=nothing.rec %s', further_arguments);
        end
        % If the 'trajectories' variable is not given at all, the workbench default .rec-file is probably used (and nowadays it is mostly set to not record).

        variables(i:i+1) = []; % Remove name,value pair from list of unprocessed variables
        
%       case 'source_point_offset_vertical'
%         % (Has already been used to produce FLY2-file, but remove from list here,
%         % because it is not an adjustable to pass to Lua sript.)
%         source_point_offset_vertical = variables{i + 1};
%         variables(i:i+1) = []; % Remove name,value pair from list of unprocessed variables
%         if source_point_offset_vertical ~= 0 && source_point_pattern ~= -1
%           error('Non-zero source_point_offset_vertical is only supported when 3D Gaussian source point distribution by fly2-file. Could be fixed by adding it as an adjustable in Lua script and making the Lua script set ion_pz_mm for any pattern where it sets ion_px_mm and ion_py_mm.');
%         end
        
      otherwise
        % Unknown entry, skip it and its value (will raise one error after loop, listing all unrecognized entries)
        i = i + 1;
    end
  end
end
if ~isempty(variables)
  disp('Got variables = ');
  disp(variables);
  error('The variable name "%s" (and possibly others) was not recognized.', variables{1});
end

%write_DLT_output = ~isempty(variables) && get_argument_from_cells(variables,'write_DLT_output');
write_DLT_output = nargout >= 3; % Let enabling DLT writing be as simple as asking for the output filename! It has to be handled anyway as it is currently written to the constant name "sim.dlt".
if write_DLT_output
  % IMPROVEMENT: make the DLT writer work, use _G (global) set from .fly2 or adjust.lua instead if the --adjustable wouldn't work
  further_arguments = [further_arguments ' --adjustable write_DLT=4']; % make the workbench Lua script write DLT header and footer
end

% command = sprintf('"%s" --nogui fly --particles=%s --recording-output=%s --retain-trajectories=0 %s.iob %s', simion_path, fly2_name, output_name, workbench_name, further_arguments);
% 2016-01-07 15:42 found in SIMION docs that loading of potentials from IOB can be skipped (saving time since we will override them) by --restore-potentials=0
command = sprintf('"%s" --nogui fly --restore-potentials=0 --particles=%s --recording-output=%s --retain-trajectories=0 %s.iob %s', simion_path, fly2_name, output_name, workbench_name, further_arguments);
% -- didn't notice any speedup, but not worse either


% DEBUG can try to handle larger potentail arrays (reflectron5 without cylinder symmetry). (With 8GB RAM up to both3.iob was OK without the next line):
% command = sprintf('"%s" --reserved-memory=3G  --nogui fly --particles=%s --recording-output=%s --retain-trajectories=0 %s.iob %s', simion_path, fly2_name, output_name, workbench_name, further_arguments);

% --- ---
% Prepare to run Simion
previous_directory = cd;
try % to return to previous directory on error
  cd(simulation_root)

  % Write to .fly2 file, using unique/random name to perhaps allow parallell work
  f = fopen(fly2_name, 'w', 'b');
  fprintf(f, '%s', fly2);
  fclose(f);

  if exist(output_name,'file') && ~isempty(strfind(output_name,'.txt')) % the .txt check is to avoid deleting file if it does not have .txt suffix
    delete(output_name)
    if exist(output_name,'file')
      cd(previous_directory);
      error(sprintf('Could not remove old sim.txt file.\nYou may need to terminate a hidden Simion process.'))
    end
  end
  if exist('adjust.lua','file')
    % Disable "adjust.lua" so that parameters are not overridden by it
    if exist('adjust -- disabled.lua','file')
      delete('adjust -- disabled.lua')
    end
    movefile('adjust.lua', 'adjust -- disabled.lua');
  end
  if ~isempty(adjustment_script)
    % Write to adjust.lua file, using unique/random name to perhaps allow parallell work
    % WARNING NOTE: this breaks parallellism, since the adjust.lua name is currently not suffixed by any counter or random string.
    f = fopen('adjust.lua', 'w', 'b');
    fprintf(f, '%s\r\n', adjustment_script);
    fclose(f);
    disp('Writing to adjust.lua.') % DEBUG
  end

  if debug_mode
    disp(command);
  end
  % Run
 	[status, output] = dos(command);
  if status ~= 0
    cd(previous_directory);
    disp(command);
    error(sprintf('Something went wrong when calling Simion. Status: %d\nOutput:\n%s', status, output));
  end
  
  % Remove (old) temporary files
  delete *.tmp
  if ~debug_mode 
    delete(fly2_name);
    if ~isempty(adjustment_script)
      delete('adjust.lua')
    end
  end

  if write_DLT_output
    DLT_path = strrep(output_name, '.txt','.dlt')
    if strcmp(DLT_path,output_name) || ~exist(DLT_path,'file')
      % no DLT output
      DLT_path = '';
    else
      % Return DLT_path so caller can copy file to a permanent storage
      DLT_path = fullfile(cd, DLT_path);
    end
  else
    DLT_path = '';
  end
  
catch e
  cd(previous_directory);
  rethrow(e);
end
cd(previous_directory);

if nargout >= 1 % Also call readsim3 to read the logged data from the simulation
  logfile_path = [simulation_root output_name];
  flys = readsim3(0, [simulation_root output_name]);
  if isstruct(flys) && ~isempty(flys)
    % readsim doesn't know the workbench name or source point (unless trajectories kept). Add it to the retained info:
    for i = 1:length(flys)
      flys(i).workbench_name = workbench_name;
      
      % Store also the source point coordinates (at least to the level known here within runsim3, assuming the same default coordinate is defined for the workbench as in in Lua script actually used) 
      flys(i).source.normal_source_point_z = normal_source_point_z; % a property of the workbench, without the possibly used offset for this simulation
      if ~isfield(flys(i).source, 'source_point_z')
        flys(i).source.source_point_z = source_point_z; % z for position is the TOF-axis (x in SIMION, z in Matlab)
      end
      if ~isfield(flys(i).source, 'source_point_y')
        flys(i).source.source_point_y = source_point_y;
      end
      if ~isfield(flys(i).source, 'source_point_vertical')
        flys(i).source.source_point_vertical = source_point_vertical;
      end
      if ~isfield(flys(i).source, 'initial_v_vertical')
        % Include jet velocity info, and translate back from the old bad name initial_vz to initial_v_vertical 
        % (z in Matlab & most of Lua-adjustable notation is the TOF-axis, not the vertical axis).
        flys(i).source.initial_v_vertical = initial_v_vertical;
      end
      
    end
    % NOTE: new Lua scripts print "Detector coordinates" which readsim3 puts in flys(end).geometry
    
    if ~debug_mode && nargout <= 2
      % If the output filename is not requested, the parsed data in flys is enough. Delete the log file.
      delete(logfile_path);
    end
  else
    % keeping logfile in case result was empty (since some error is likely to need debugging)
    flys = []; % use empty array instead of NaN (if that is still returned by readsim on error)
  end
end
