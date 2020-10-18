% Simulate and show VMI images.

constants;
first_job = 1;

% A list of settings to iterate over, so the program can run a sequence of simulations without user interaction
names_workbenches_and_potentials = { 
'' 'example', [-3000,-2095,   0]; % electron VMI. Ratio 0.698, optimized manually with nominal source point.
% '' 'example',[-3000,-2140,   0]; % electron VMI. Ratio 0.7133.
% '' 'example',[-4600,-3308,   0]; % electron VMI tuned for increased energies
% '' 'example',[4000,3400,0]; % ions, mass spectrum
};

potential_offset_components = { % a basis set for trying deviations, e.g. of +-50 V.
};
potential_offset_values = []; % don't do any voltage offsets
% potential_offset_values = [-50 0 50]; % original and 50 V in either direction. Needs 3^size(potential_offset_components,1) evaluations if fully combinatorial.

if isempty(potential_offset_components)
  potential_offset_components_sequential = [];
  potential_offset_components_combinatorial = [];
  potential_offset_names_sequential = {};
  potential_offset_names_combinatorial = {};
else
  potential_offset_components_sequential    = cell2mat(potential_offset_components(~cell2mat(potential_offset_components(:,1)),2));
  potential_offset_components_combinatorial = cell2mat(potential_offset_components( cell2mat(potential_offset_components(:,1)),2));
  potential_offset_names_sequential    = potential_offset_components(~cell2mat(potential_offset_components(:,1)),3);
  potential_offset_names_combinatorial = potential_offset_components( cell2mat(potential_offset_components(:,1)),3);
end
if ~isempty(potential_offset_components_combinatorial)
  %error('potential_offset_combinatorial=true is not implemented, may be possible if only a few components selected by true in potential_offset_components{:,1}.');
  combinatorial_count = length(potential_offset_values) ^ size(potential_offset_components_combinatorial,1);
  if length(potential_offset_values) > 3
    error('Combinatorial voltages with more than %d values is not supported.', length(potential_offset_values));
  end
  fprintf('COMBINATORIAL (%d values)^(%d directions) : %d evaluations\n', length(potential_offset_values), size(potential_offset_components_combinatorial,1), combinatorial_count);
else
  combinatorial_count = 1;
end
if ~isempty(potential_offset_components_sequential)
  fprintf('Sequential potential_offset_components of %d values in in %d directions: %d evaluations\n', length(potential_offset_values), size(potential_offset_components_sequential,1), size(potential_offset_components_sequential,1)*length(potential_offset_values));
end
if ~isempty(potential_offset_components_combinatorial) || ~isempty(potential_offset_components_sequential)
  if isempty(potential_offset_values)
    error('Can''t run combinatorial or sequental voltage offsets when the list of offsets is empty.');
  end
elseif ~isempty(potential_offset_values)
  error('If the list of voltage offsets components is empty, also the potential_offset_values should be empty.');
end
reflectron_log_inputs = {}; % workbench_name, job_name, offsets, potential_offset_str, potential_offset_name, potential_offset_amount, potentials_orig, potentials_offset, potentials, variables
reflectron_log_results = []; % mean(H,V,t), std(H,V,t), [H.linear, H.cubic, V.linear, magnification_radial], [t/std(t), mass overlap[%], dT/dz0(linear term)]

% Command for importing voltage in SIMION GUI with };_G.U[0]=0;U=_G.U; adjustable DoNotSetPotentials=-1; _G.FLY2_IONS_PER_GROUP=32
%  sprintf('_G.U = {%s};', num2str(potentials, '%.8g, '))

%% Source offsets
offset_list = [ % Offsets along: optical, vertical, spectrometer(TOF) axes [mm], and upwards velocity (-jet) [km/s]
  0, 0 , NaN, NaN % centered
%   0, -0.5, NaN, NaN % (small) vertical offset (down, in same direction as jet velocity)
%   0, 1, NaN, NaN % vertical offset
  % 1, 0, NaN, NaN % small horizontal offset
  % 0, 2, NaN, NaN % (large) vertical offset
%   2, 0, NaN, NaN % horizontal offset
%   4, 0, NaN, NaN % (large) horizontal offset
  
  % To look for asymmetry due to source point offset along optical axis (y in SIMION but shown horizontal as in lab):
%   8, 0, NaN, NaN %[mm] huge offset (near extractor radius 10 mm)
  % To run with additional z-offset (vertical) to see if any combined effect (e.g. together with camera tilt) that makes apparent C-shape deviations
% 	8, 2, NaN, NaN %[mm] to test combined effect of offsets in both directions. Seems to show that this is the case.

  5, 1, NaN, NaN % A smaller offset variant

% Controlling also initial_v_vertical by this list:
%   0, 0 , NaN, 0 % centered, no jet velocity
 0, 0 , NaN, -0.3 % MAIN CHOICE: centered, 0.3 km/s jet downwards. 
 % As long as there is some vertical source point spread, this setting gives all the interesting info 
 % (similar spatial mapping coefficients as if run with vertical offset or withou jet, but with the 
 % info about jet-mass-separation). Speed of sound in air is about 0.34, in xenon about 0.17)
];

%%%%%%%%%%%%%%%%%%%%%%%%%
show_starts = false; % normal simulation
% show_starts = true;  % to visualize starting points in same coordinate system

for job = 1:size(names_workbenches_and_potentials,1)
  job_name = names_workbenches_and_potentials{job,1};
  job_name = strrep(strrep(job_name, '-',''), '+',''); % OPTION to hide any +- suffixes that just indicate how nice the setting is
  workbench_name = names_workbenches_and_potentials{job,2};
  potentials_orig = names_workbenches_and_potentials{job,3};
  detectorB_linecoeff = [];
  
  % Which indices in potentials_orig correspond to the experimentally main controlled for VM?
  switch workbench_name
    case 'example'
      pot_index_E12 = 3;
      pot_index_E14 = 1;
      pot_index_E15 = 2;
      pot_index_E18 = 3;
      pot_index_deflectors = [];
      detector_position = 283; %[mm] detectorA_x = 947 -- A-MCP front surface position (depends on the geometry).
    case 'ecyl8e'
      pot_index_E12 = 2;
      pot_index_E14 = 4;
      pot_index_E15 = 5;
      pot_index_E18 = 8;
      pot_index_deflectors = [];
      detector_position = 397; %[mm] detectorA_x = 397 -- A-MCP front surface position (depends on the geometry).
    case {'both8eu', 'both8eu_B2'}
      pot_index_E12 = 12;
      pot_index_E14 = 14;
      pot_index_E15 = 15;
      pot_index_E18 = 18; 
      if strcmp(workbench_name, 'both8eu_B2')
        if length(potentials_orig) ~= 19
          error('both8eu_B2 requires 19 potentials.')
        end
      else
        if length(potentials_orig) ~= 18
          error('both8eu requires 18 potentials.')
        end
      end
      pot_index_deflectors = [7 9 10 11];
      detector_position = 947; %[mm] detectorA_x = 947 -- A-MCP front surface position (depends on the geometry).
      % detectorB_linecoeff is Lua's detectorB_linecoeff = {k, -1, -d, MCP_centre_y}:
      % k=-1/tan(det.angle=6deg), d = k*MCP_centre_x - 1*MCP_centre_y; using MCP imported as separate PA
      detectorB_linecoeff = [-9.5144, -1.0000, 3805.5447, 133.00];
    otherwise
      error('TODO implement');
  end
  
  mass_in_u = electron_mass / atomic_mass_unit; %[u]
  %% Array of kinetic energies and masses (&charges) to simulate for.
  %original_energies = [0 1 2 3 4]; %[eV] having an exact zero gives an extremely strong intensity, making weak tails visible
  %original_energies = [k_Boltzmann*(273+19)/eV 1 2 3 4]; %[eV] The central "spot" is quite big compared to experiment, suggesting jet is cooled below room temperature
  %original_energies = [k_Boltzmann*(150)/eV 1 2 3 4]; %[eV] The central "spot" is quite small in experiment, suggesting jet is cooled below room temperature. Used for UV experiment settings when max 10 eV
  original_energies = [k_Boltzmann*(70)/eV 1 2 3 4]; %[eV] The central "spot" is quite small in experiment, suggesting jet is cooled below room temperature
  % original_energies = [0.1]; %[eV]
%   original_energies = [4]; %[eV]
  %original_energies = [k_Boltzmann*(150)/eV, [0.7 2], [0.7 2]+1.7]; %[eV] The central "spot" is quite big compared to experiment, suggesting jet is cooled below room temperature
  %original_energies = [4.66 6.23 6.81 7.80 9.37 10.94 12.51 14.08 15.65 17.22 18.79]; %[eV] He and H2O by long harmonics of 790 nm, to order 20
%   original_energies = [15]; %[eV] just for testing the peak shape after Abel inversion, effects due to extent of source volume (non-centered source trajetories get higher lensing and smaller radius and lower apparent kinetic energy: asymmetric peak shape)
  original_energies = [0.3, 1, 2, 4, 6, 10]; %[eV] for VMI calibration 2020-09-15
%   original_energies = [0.3, 1, 3, 6, 10, 15, 20, 25, 30]; %[eV] for VMI calibration 2020-09-15
  
  if potentials_orig(end) < -1000  % seems like VMI in ion mode
    % To simulate cations
    mass_in_u = 131.2; %[u] standard atomic weight of xenon
    mass_in_u = 28; %[u] N2
    %original_energies = [k_Boltzmann*(150)/eV, k_Boltzmann*(273+20)/eV, 2]; %[eV]
    original_energies = [k_Boltzmann*(80)/eV, k_Boltzmann*(273+20)/eV, 2]; %[eV]
%     original_energies = [k_Boltzmann*(80)/eV, k_Boltzmann*(273+20)/eV, 1, 2, 4, 6]; %[eV] for VMI calibration 2020-09-15
  elseif potentials_orig(1) > 1000 
    % To simulate cations (water in Zuerich_guess.iob)
    mass_in_u = 18; %[u] water
    original_energies = [k_Boltzmann*(80)/eV, k_Boltzmann*(273+20)/eV]; %[eV]
  end
%   disp('Using a list of ions.'); % OPTION to just simulate ions:
%   mass_in_u = [18 18 40 40]; %[u] water and argon cations
%   mass_in_u = [2 2 131 131]; %[u] dihydrogen(included in spatial dependence plot) and xenon, to get something with large mass difference
% %   mass_in_u = [1 1 131 131]; %[u] hydrogen atom (excluded from spatial dependence plot) and xenon, to get something with large mass difference
% %   mass_in_u = [40 40]; %[u] argon cations
% %   mass_in_u = [18 18]; %[u] water cations
% %   original_energies = repmat([k_Boltzmann*(30)/eV, k_Boltzmann*(273+20)/eV], 1, ceil(length(mass_in_u)/2)); %[eV]
%   original_energies = repmat(k_Boltzmann*(80)/eV, 1, length(mass_in_u)); %[eV] the same (low) energy for all ions
% %   original_energies = repmat(1, 1, length(mass_in_u)); %[eV] the same (high) energy for all ions, just to see the effect of such spread (as if fragment from molecule with large KER)
%   
%   %mass_in_u=[1 3000]; original_energies=[1 1]*k_Boltzmann*(80)/eV; % 80 K for hydrogen as well as heavy ion, to try to resolve spots for mass-jet-shift (quite small effect, smaller than assumed vertical source size)
%   %mass_in_u=[1000]; original_energies=[1*k_Boltzmann*(80)/eV]; % a single heavy ion to compare with next (see that 999 overlaps and outer halo is dihydrogen)
%   mass_in_u=[2 999 1000]; original_energies=[1, [1 1]*k_Boltzmann*(80)/eV]; % 1eV for dihydrogen, then 80 K for the heavy ion pair to show some mass resolution too. Good for reflectron overview.
% %   mass_in_u=[1 18 40]; original_energies=[1, [1 1]*k_Boltzmann*(80)/eV]; % 1eV hydrigen, then 80K water and Ar

%   mass_in_u=[18]; original_energies=[k_Boltzmann*(20)/eV]; % just water, and quite cold
%   mass_in_u=[18]; original_energies=[k_Boltzmann*(273)/eV]; % just water, room temperature (ignoring any prefactor like 3/2 or 2/3)

  
  %% Source size (source_point_distribution)
  % Use a pseudorandom (reproducible) 3D gaussian with these standard deviations
  % along the SIMION X,Y,Z coordinates. If std_transverse_other is omitted it defaults tostd_transverse.
  % If an empty array is given, the default is used. Converting from FWHM [mm] along (TOF axis, optical axis, vertical (jet) axis).
  % Default: currently [0.06 1 0.06] 3D gaussian.
  source_point_pattern = [0.2 8 0.2] / (2*sqrt(2*log(2))); % 8 mm long and not very narrow
  
  for offset_index = 1:size(offset_list,1)
    source_point_offset_y = offset_list(offset_index, 1);
    source_point_offset_vertical = offset_list(offset_index, 2);
  
  combinatorial_indices = ones(1,size(potential_offset_components_combinatorial,1));
  for dV_combinatorial_index = 1:combinatorial_count
  for dV_seq_index = 1:max([1 size(potential_offset_components_sequential,1)])
  for dV_index = 1:iif(isempty(potential_offset_components_sequential), 1, length(potential_offset_values))
    potential_offset_name = '';
    potential_offset_component = [];
    potential_offset_amount = 0;
    potential_offset = 0;
    potential_offset_str = '';
    if ~isempty(potential_offset_components_combinatorial)
      % Implementing combinatorial offsets somehow...
      %potential_offset = 0;
      potential_offset = potential_offset_values(combinatorial_indices) * potential_offset_components_combinatorial;
      potential_offset_str = '';
      for cm = 1:length(combinatorial_indices)
        % Use compact notation: a single letter, so the position encodes voltage 
        % and for one of the values (the second) the single-letter name of the component is shown.
        if combinatorial_indices(cm) == 1
          % For first value (usually negative), use a low ASCII value
          potential_offset_str = [potential_offset_str '#'];
        elseif combinatorial_indices(cm) == 2
          % Second value (which is normally positive if two, or 0 if three)
          potential_offset_str = [potential_offset_str potential_offset_names_combinatorial{cm}];
        else % third value, use a high ASCII value to sort last
          potential_offset_str = [potential_offset_str '¶'];
        end
      end
      % Include one voltage (assuming all are +- (or 0) by the same amount
      potential_offset_str = sprintf('%s(%.0fV)', potential_offset_str, max(abs(potential_offset_amount)));
    end
    if ~isempty(potential_offset_components_sequential) % Sequence of independentent potential offsets
      potential_offset_name = potential_offset_names_sequential{dV_seq_index};
      potential_offset_component = potential_offset_components_sequential(dV_seq_index,:);
      potential_offset_amount = potential_offset_values(dV_index);
      if isempty(potential_offset)
        potential_offset = potential_offset_component * potential_offset_amount;
      else
        potential_offset = potential_offset + potential_offset_component * potential_offset_amount;
      end
      if isempty(potential_offset_str)
        potential_offset_str = sprintf('%s%+.0fV', potential_offset_name, potential_offset_amount);
      else
        potential_offset_str = sprintf('%s,%s%+.0fV', potential_offset_str, potential_offset_name, potential_offset_amount);
      end
    end
    potentials = potentials_orig + potential_offset;
    
    if strcmp(workbench_name, 'both8eu_B2')
      % Show the experimentally often tweaked E10bottom (index 7) and E9push
      %pot_str = sprintf('%s E10b,9;12,14,15=%.0f,%.0f;%.0f,%.0f,%.0fV', workbench_name, ...
      pot_str = sprintf('%.0f,%.0f;%.0f,%.0f,%.0fV', ...
          potentials(7), potentials(9), potentials(pot_index_E12), potentials(pot_index_E14), potentials(pot_index_E15));
    elseif potentials(pot_index_E18) == 0
      pot_str = sprintf('E12,14,15=%.0f,%.0f,%.0fV', potentials(pot_index_E12), potentials(pot_index_E14), potentials(pot_index_E15));
    else
      pot_str = sprintf('E12,14,15,18=%.0f,%.0f,%.0f,%.0fV', potentials(pot_index_E12), potentials(pot_index_E14), potentials(pot_index_E15), potentials(pot_index_E18));
    end
    if ~isempty(job_name)
      % If a job_name is provided, assume it is unique so that job_name + potential_offset_str is enough to show (not any of the resulting potential values)
      pot_str_w_offset = potential_offset_str;
    else
      % If no job_name, include both offset and some of the resulting potentials
      pot_str_w_offset = [potential_offset_str ' ' pot_str];
    end
    if ~strcmp(workbench_name, 'both8eu_B2')
      % include name only if not the usual full model
      pot_str = sprintf('%s %s', workbench_name, pot_str);
      pot_str_w_offset = sprintf('%s %s', workbench_name, pot_str_w_offset);
    end
    if isempty(potential_offset_str)
      % When no potential offset is used (normal case)
      pot_str_w_offset = pot_str;
    end
    
    % from optimize_potentials.m for 'spheres' option:
    multipole = 0; % isotropic
    %multipole = 1; % dipole = |first order Legendre polynomial|^2
%     multipole = 3; % |third order Legendre polynomial|^2. A bit more aligned along polarization axis
    % multipole = 4; % |fourth order Legendre polynomial|^2
%     multipole = 7; % |seventh order Legendre polynomial|^2. Used for PECD campaign UV & IR

    if length(original_energies) <= 2
      n = 20000; % high statistics for few rings
      n = 50000; % high statistics for a single ring when small pixel size
    else
      n = 10000; % (still quite a lot)
    end
%     n = 20000; % high statistics despite many energies
    n = 5000; % quite OK for testing, not for publication
%     n = 1000; % for quick testing, not for publication
%     n = 500; % for quick testing, not for publication
    
    if show_starts
      if length(original_energies) == 1
        n = 4000; % for quick testing, not for publication
      else
        n = 1000; % for quick testing, not for publication
      end
    end
    
    % Isotropic or dipole-transition like (density proportional to z^2), random directions, no grouping for VMI resolution calculation:
    variables = {'_number_of_directions',2, 'multipole_order',multipole, 'particles_per_direction',ceil(n/2), '_isotropic',1, 'trajectories',true};
    if multipole == 0
      % Fully isotropic (not alignment along polarization axis)
      variables = {'_number_of_directions',1, 'multipole_order',multipole, 'particles_per_direction',ceil(n/2), '_isotropic',1, 'trajectories',true};
    end

%     variables{1,end+1}='debug'; variables{1,end+1}=true; % -- retain the job files afterwards instead of deleting them when runsim3 had read it

    source_point_offset_z = offset_list(offset_index, 3); % to optionally specify via offset_list and its loop
    if ~isempty(source_point_offset_z) && isfinite(source_point_offset_z)
      variables{1,end+1}='source_point_offset_z'; variables{1,end+1} = source_point_offset_z; % [mm] along spectrometer (TOF) axis
    end
    initial_v_vertical = offset_list(offset_index, 4); % to optionally specify via offset_list and its loop
    if ~isempty(initial_v_vertical) && isfinite(initial_v_vertical)
      variables{1,end+1}='initial_v_vertical'; variables{1,end+1} = initial_v_vertical; % [km/s] molecular beam velocity approx for helium (liquid water 1.5 but vapor probably lower). Insignificant velocity offset for electrons.
    else
      % Old way of setting initial_v_vertical was to edit here, but can use offset_list instead now
      % variables{1,end+1}='initial_v_vertical'; variables{1,end+1} = -0.3; disp('With jet velocity.'); % constant offset in intial velocities [mm/us]=[km/s] e.g. to represent speed of molecular beam -- only added to Lua program for some geometries!
    end
    variables{1,end+1}='source_point_offset_y'; variables{1,end+1} = source_point_offset_y;
    variables{1,end+1}='source_point_offset_vertical'; variables{1,end+1} = source_point_offset_vertical;
    
    % Fully customized potentials, given as a row array following the SIMION indexing of electrodes.
    variables{1,end+1}='potentials'; variables{1,end+1}=potentials;
    
    fprintf('Running %d*%d with offset (%.4g, %.4g, %.4g) mm (or z0=%.4g) for %s. %s\n', n, length(original_energies), ...
      source_point_offset_y, source_point_offset_vertical, source_point_offset_z, ...
      get_argument_from_cells(variables, 'source_point_z', 'num', NaN), pot_str, ...
      mat2str(source_point_pattern*(2*sqrt(2*log(2)))));
    
    disp(strcat('_G.U={', num2str(potentials,'%.0f,'), '};_G.U[0]=0;U=_G.U; adjustable DoNotSetPotentials=-1;'))
    
    % Run the simulation (can take many minutes when more than 1000 particles)
    [flys,logfile_path] = runsim3(workbench_name, mass_in_u, original_energies, source_point_pattern, variables);
    
    
    %% Visualize the simulated results
    % This cell can be run also just after loading a .mat file with saved simulation trajectory data.
    % NOTE if loading old .mat file without all variables, use potentials=get_argument_from_cells(variables, 'potentials','num')
    
    figure(3); clf;
    set(3, 'Position', [1055 59 748 698]);
    pair_counter = 0;
    mass = collect_field({flys.source}, 'mass', 2);
    miss_fraction = sum(collect_field(flys, 'misses', 2)) / sum(collect_field({flys.source}, 'ion_count', 2));
    
    image_from_detector_B = false;
    if all(mass_in_u >= 1)
      % Runnig ion instead of electron
      particle_indices = find(mass >= 1);
      if potentials_orig(end) < -1000  % seems like VMI in ion mode
        image_from_detector_B = false;
      else % Assume that reflectron detector is used
        image_from_detector_B = true;
      end
    else % Electrons as expected for VMI
      particle_indices = find(mass < 1);
    end
    if isempty(particle_indices)
      error('No particles selected for display!')
      %continue; % nothing to show
    end
    colours = 'ckbgr';
    colours = 'kybgr';
    colours = 'krmgb';
    if length(particle_indices) == 2
      colours = 'gb';
    end
    
    detectorA_centre = [0, 0]; % [mm] optical axis, vertical axis
    detectorB_centre = [NaN, NaN]; % [mm] optical axis, vertical axis
    % Normal case, VMI-side detector ("A")
    if isfield(flys(end),'geometry') && ~isempty(flys(end).geometry)
      detectorA_centre = flys(end).geometry.detector_A_coordinates(2:3); % [mm] along optical axis and vertical axis
      if detector_position ~= flys(end).geometry.detector_A_coordinates(1)
        error('Conflicting values for detector_position along TOF axis: %.2f from Lua and %.2f defined here', flys(end).geometry.detector_A_coordinates(1), detector_position);
      end
    end
    if image_from_detector_B
      % New option of using image from ion-side reflectron detector ("B")
      if ~isfield(flys(end),'geometry') || isempty(flys(end).geometry)
        error('Position information for detector_B position is not available from the read simulation log. It is required for using detector "B" (reflectron).');
      end
      if isempty(detectorB_linecoeff)
        % Simply get the values from Lua
        detectorB_linecoeff = flys(end).geometry.detectorB_linecoeff;
      elseif any(detectorB_linecoeff ~= flys(end).geometry.detectorB_linecoeff)
        % Values already given, but conflict with Lua
        error('Conflicting values for detector_B position and orientation: %s from Lua and %s defined here', mat2str(flys(end).geometry.detectorB_linecoeff), mat2str(detectorB_linecoeff));
      end
      % Currently always using flys(end).geometry.detectorB_xSIMION_zMatlab and flys(end).geometry.detector_A_coordinates(3) from Lua, not defining them here too.
      detectorB_centre = [detectorB_linecoeff(4), flys(end).geometry.detector_A_coordinates(3)]; % [mm] along optical axis and vertical axis
      detector_name = 'Reflectron';
    else
      detector_name = 'VMI';
    end
    mean_energy = mean(original_energies);
    
    if flys(particle_indices(1)).source.ion_count < 1000 && ~show_starts
      image_bin_width = 2; % [mm]
    elseif flys(particle_indices(1)).source.ion_count <= 10000
      image_bin_width = 0.4; % [mm]
    else
      image_bin_width = 0.2; % [mm]
      image_bin_width = 0.25; % [mm]
      if n > 30000
        image_bin_width = 0.082; %[mm/px] using experimentally determined value for camera of 2019-10 water campaign
      end
    end

    halfwidth = 36; %[mm] show up to circa +-36 mm x and y 
    image_centre_index = 1 + ceil(halfwidth / image_bin_width);
    image_bin_max = 2*image_centre_index-1; % the indices from 1 to image_centre_index-1 represent negative coordinates, image_centre_index+1 to 2*2*image_centre_index-1 represent positive
    image_histogram = zeros(image_bin_max,image_bin_max);
    one = uint32(1);
    
    % Define variables that accumualate results for all the selected runs (masses & charges & energies)
    legend_strings = {};
    energies = [];
    h_for_legend = [];
    colours=fliplr(colours); particle_indices = fliplr(particle_indices); % to draw lowest energy dots last (on top)
    radii = [];
    coords_bad_all = [];
    max_shown_points = 4000;
    angular_halfwidth_selection_for_magnification = NaN;
    angular_halfwidth_selection_for_VMI = NaN;
    plane_radii_and_energies = []; % for energy calibration, using only elevation angles near 90 degrees
    magnifications_and_start_radii = NaN(0,2); % for spatial mode magnification (using only trajectories with minimal transverse velocity)
    start_and_final_coords_where_nearly_axial_velocity = NaN(0,4); % for axis resolved spatial dependency plot (using only trajectories with minimal transverse velocity)
    TOF_and_start_z0_all = NaN(0,3); % time of flight [ns], initial position along TOF axis [mm], mass [u].
    
    for i = 1:length(particle_indices) % for electrons (in energy-reversed order)
      f = flys(particle_indices(i)); % results for this energy
      
      % NOTE: X and Z are swapped from SIMION coordinates to the readout names in readsim3.
      % In SIMION: x is the TOF axis but it is called trajectories.z here,
      % In SIMION: z is the vertical axis (polarization) but it is called trajectories.x here. Shown vertically in diagrams. Positive means up in reality (although one SIMION GUI diagram is drawn with it pointing down).
      % In SIMION: y is the optical axis, it is called trajectories.y. Shown horizontally in diagrams here.
      % Assume that (as still in all models) the "A" VMI-detector is centered on the nominal source point.
      coords_start = [f.trajectories.y0, f.trajectories.x0] - detectorA_centre; %[mm] source coordinates. First horizontal (optical axis), then vertical.
      axial_start = f.trajectories.z0 - f.source.normal_source_point_z; %[mm] along the TOF-axis
      t = f.trajectories.t; %[ns] time of flight
      if ~image_from_detector_B
        % Normal case, VMI-side detector ("A")
        % NOTE: when plotting VMI image with optical axis horizontal, it means y is horizontal (optical axis) and x vertical...
        coords = [f.trajectories.y, f.trajectories.x] - detectorA_centre; % detector coordinates. First horizontal (optical axis), then vertical.
        bad = abs(f.trajectories.z - detector_position) > 3; % reject particles not impacting the detector (e.g. hitting the extractor electrode)
      else
        % New option of using image from ion-side reflectron detector ("B")
        % detectorB_linecoeff is [k, -1, -d, MCP_centre_y]
        %   k = -1/tan(det.angle=6deg),
        %   d = k*MCP_centre_x - 1*MCP_centre_y; %(MCP_centre_x=nominal_z=386 as the MCP import Xwb)
        % Detector surface is a line when projected onto (z,y)-plane. Check that deviation from that line is smaller than z_margin
        % reject particles not impacting the detector (e.g. hitting the extractor electrode)
        bad = abs( detectorB_linecoeff(1)*f.trajectories.z  + detectorB_linecoeff(2)*f.trajectories.y + detectorB_linecoeff(3) ) > 1.5; %[mm], not using any max_r. Tolerance as in show_mass_resolution.m. A test run with 1000 particles shows a tolerance of 0.02 mm would usually be enough, but a few millimetres gives no artefacts.
        
        nominal_y = detectorB_linecoeff(4);
        if detectorB_linecoeff(1) ~= 0 % normal case, the angle is nonzero (detector surface normal is not parallel to the z-axis)
          nominal_z = (-detectorB_linecoeff(2)*nominal_y - detectorB_linecoeff(3)) / detectorB_linecoeff(1);
        else %special case when detector normal is along z axis (in this case the roles of y and z axes become somewhat reversed to what is assumed mostly)
          nominal_z = detectorB_linecoeff(4);
          nominal_y = (-detectorB_linecoeff(1)*nominal_z - detectorB_linecoeff(3)) / detectorB_linecoeff(2);
        end
        % Project onto the detector surface (it is rotated 90 degrees wrt. the normal vector in detectorB_linecoeff(1:2))
        % Long form:
        %   vector_projected_onto_normal         = detectorB_linecoeff(1:2)' * detectorB_linecoeff(1:2)*[(f.trajectories.z'-nominal_z);(f.trajectories.y'-nominal_y)]/sum(detectorB_linecoeff(1:2).^2);
        %   vector_projected_onto_detector_plane = [(f.trajectories.z'-nominal_z);(f.trajectories.y'-nominal_y)] - vector_projected_onto_normal; % get the remaining vector in detector plane
        %   coords = [transpose([detectorB_linecoeff(2) -detectorB_linecoeff(1)] * vector_projected_onto_detector_plane / norm(detectorB_linecoeff(1:2))), f.trajectories.x-detectorB_centre(2)];
        % Short form, giving same result:
        coords = [(detectorB_linecoeff(2)*(f.trajectories.z-nominal_z) - detectorB_linecoeff(1)*(f.trajectories.y-nominal_y)) / norm(detectorB_linecoeff(1:2)), f.trajectories.x-detectorB_centre(2)];
        % coords(:,1) is horizontal coordinate in detector's (rotated) coordinate system (mixture of optical axis and TOF axis).
        % coords(:,2) is vertical axis, unchanged.
        % flys(end).geometry.detectorB_xSIMION_zMatlab; % not sure this is correct, but not needed here
      end
      
      coords_bad = coords(bad,:);
      if show_starts
        coords = coords_start; % DEBUG to show source coordinates instead of result of trajectory
        bad = false & bad;
      end
      
      % Reject particles not impacting the detector (e.g. hitting the extractor electrode)
      coords(bad,:) = []; 
      coords_start(bad,:) = [];
      axial_start(bad) = [];
      t(bad) = [];
      coords_bad_all = [coords_bad_all; coords_bad];
      if mass(particle_indices(i)) >= max(mass) - 1.1
        % Only log the TOF for the heaviest mass (including heaviest-1u if present) (possibly multiple energies), to get single TOF reference point.
        TOF_and_start_z0_all = [TOF_and_start_z0_all; t, axial_start, repmat(f.source.mass, length(t), 1)]; % [ns], [mm] 
      end
      
      if size(coords,1) > max_shown_points
        % Plot only the first 4000 poins, since it becomes hard to see intensity variations with too many points
        h = plot(coords(1:max_shown_points,1), coords(1:max_shown_points,2), ['.' colours(1+mod(i-1,length(colours)))], 'MarkerSize',4.5);
      elseif size(coords,2) >= 2 && size(coords,1) >= 1
        h = plot(coords(:,1), coords(:,2), ['.' colours(1+mod(i-1,length(colours)))], 'MarkerSize',4.5);
      else % no points
        h = 0;
      end
      h_for_legend = [h_for_legend; h(1)];
      hold on;
      %plot(coords_bad(:,1), coords_bad(:,2), ['x' colours(i)], 'MarkerSize',2);
      %plot(coords_bad(:,1), coords_bad(:,2), ['x' 'c'], 'MarkerSize',2);
      legend_strings =  [legend_strings; sprintf('%seV, %d misses', texformat_SI(f.source.energy,2), f.misses)];  
      %energies = [energies, f.source.energy]; % Seems OK to use f.source.energy, has four decimals in eV. 
      energies = [energies, mean([f.source.energy f.trajectories.energies(1)])]; % For random/thermal energy (not a round value in eV) could possibly combine the per-trajectory value (more decimals) 
      %energy_deviation = ([f.source.energy f.trajectories.energies(2)] - original_energies(find_nearest(original_energies, f.source.energy, 1E-3))) / 1E-6 % DEBUG
      
      
      bins = round(coords / image_bin_width) + image_centre_index;
      % Accumulate into 2D histogram (image)
      valid = all(bins >= 1 & bins <= image_bin_max, 2); % skip particles outside selected image size
      bins(~valid,:) = 1; % use index [1;1] for invalid ones to avoid error, their increment will be zeroed anyway so they don't affect image
      increment = ones(size(bins,1), 1);
      increment(~valid) = 0; % don't count invalid hits into the histogram
      % NOTE: the coordinates need to be tranposed (or the bins-array-columns swapped) to get vertical coordinate (coords(:,2)) in first matrix dimension
      image_histogram = image_histogram + transpose(accumarray(bins, increment, image_bin_max*[1 1])); % make 2D-histogram without explicit for-loop

      % Accumulate list of radii
      radii_here = hypot(coords(:,1),coords(:,2));
      radii = [radii; radii_here];
      
      elevations = f.trajectories.elevation(~bad); 
      %NOTE: the elevations include the perturbation by initial_v_vertical (at least in the 
      % usually simulated istropic mode), thus filtering by elevations is screwed up
      % for low-energy ions (where a velocity of 0.3 km/s is much compared to thermal energy).
      % TODO: include initial_v_vertical in f.trajectories in readsim3.m, and either there or here 
      % subtract initial_v_vertical to define an unpergurbed_elevation to keep among trajectory variables.
      
      angular_halfwidth_selection_for_VMI = 10; %[deg]
      % For energy calibration
      % only when v_z =approx= 0 (elevation approx 90 deg.) to get correct ring without Abel inversion
      which = find(abs(elevations-90) < angular_halfwidth_selection_for_VMI); % filter within +-angular_halfwidth_selection_for_VMI deg of plane
      % This first version used nominal kinetic energy, which is slightly wrong when angle is not exactly 90 degrees
      % (as some of the energy is then out of plane and then SHOULD not appear in the VMI image plane).
      % The resulting calibration coefficient goes higher than the ideal when the filtering width for elevation is increased (tried 3, 10, 30 degrees half range)
      %  plane_radii_and_energies = [plane_radii_and_energies; radii_here(which) repmat(energies(end),length(which),1)]; % mm,eV
      % Use only the v_xy (project to velocity in plane by using *sind(elevation)) and convert to kinetic energy.
      % This gives a calibration result almost independent of elevantion filtering width (angular_halfwidth_selection_for_VMI):
      plane_radii_and_energies = [plane_radii_and_energies; radii_here(which) repmat(energies(end),length(which),1) .* sind(elevations(which)).^2]; % mm,eV
      
      % don't use the highest energies when trying to define a spatial magnification (ideally just zero energy)
      % The 1.01 factor and =1E-5 are just to give a margin for rounding errors (e.g. finite number of digits in log file)
      if f.source.energy <= 2 && f.source.energy <= 1.01*mean_energy+1E-5 % [eV] 
        % For spatial mode magnification, using only elevation angles near 0 and 180 degrees
        angular_halfwidth_selection_for_magnification = 5; %[deg]
        if f.source.energy < 0.03 && f.source.mass > 1 % [eV], [u] 
          % Far sub-thermal energy and mass as for ion heavier than proton, 
          % then energy is negligible for estimating magnification and we don't need to filter by eleveation.
          %angular_halfwidth_selection_for_magnification = inf;
          angular_halfwidth_selection_for_magnification = 60; %[deg] some arbitrary large angle to still filter a bit
        elseif n <= 5000
          %angular_halfwidth_selection_for_magnification = 10; disp('Wider angular tolerance for spatial magnification calculation');
          angular_halfwidth_selection_for_magnification = 20; disp('Wider angular tolerance for spatial magnification calculation');
        end
        if get_argument_from_cells(variables, 'initial_v_vertical', 'num', 0) == 0
          % When not perturbed by initial_v_vertical. (TODO: make readsim3 provide unperturbed_elevantion)
          % Filter within a few degrees from elevation angles near 0 and 180 degrees (5 degrees gives no points if only 500 particles run, need about 5000)
          which = find(abs(elevations-90) >= (90-angular_halfwidth_selection_for_magnification));
        else
          which = ':'; % no angular filtering in case a nonzero initial_v_vertical is used
          angular_halfwidth_selection_for_magnification = inf;
        end
        radii_here = hypot(coords_start(which,1),coords_start(which,2));
        % Project detected position onto the axis defined by source point, to not require spread only along one axis,
        % then divide one extra time by initial radius (i.e. squared) to get dimensionless magnificiation factor.
        % Incorrect before 2020-04-15: magnifications = sum(coords(which) .* coords_start(which),2) ./ radii_here.^2; % not using second dimension, i.e. just optical axis, which typically probably was similar to correct.
        magnifications = sum(coords(which,:) .* coords_start(which,:),2) ./ radii_here.^2;
        magnifications_and_start_radii = [magnifications_and_start_radii; magnifications, radii_here]; % [dimensionless], [mm]
        start_and_final_coords_where_nearly_axial_velocity = [start_and_final_coords_where_nearly_axial_velocity; coords_start(which,:), coords(which,:)]; %[mm] horizontal&vertical start, horizontal&vertical detected
      end
    end
    
    % plot(coords_bad_all(:,1), coords_bad_all(:,2), ['x' 'c'], 'MarkerSize',2);
    h_bad = plot(coords_bad_all(1:30:end,1), coords_bad_all(1:30:end,2), 'xy', 'MarkerSize',4);
    if show_starts
      % Draw a ring for the extractor electrode's inner diameter
      h = rectangle(gca, 'Position',[-10 -10 20 20 ], 'Curvature',[1 1], 'EdgeColor','k', 'LineStyle',':');
    end
    
    axis([-halfwidth halfwidth -halfwidth halfwidth]);
    daspect([1 1 1]);
    grid on
    if ~isempty(coords_bad_all)
      h_for_legend = [h_for_legend; h_bad];
      legend_strings{end+1,1} = 'Not on detector';
    end
    legend(legend_strings)
    if f.source.mass < 1 && f.source.charge == -1
      particle = 'electron';
    else
      %particle = sprintf('%.0f u ion^{%+.0f}', f.source.mass, f.source.charge);
      particle = sprintf('%s u ion^{%s}', strrep(mat2str(unique(mass)), ' ', ','), sprintf('%+.0f',unique(collect_field({flys.source}, 'charge', 2))));
    end
    if size(f.rays,1) == 1
      % (Guess that some kind of isotropic was used, the variables are not available from log file)
      directions = 'isotropic distribution';
    elseif size(f.rays,1) == 2
      % (Guess that _isotropic=1 was used with _number_of_direction=2 to get dipole distribution,
      %  the variables are not available from log file)
      switch multipole
        case 1
          directions = 'dipole(P1) distribution';
        case 3
          directions = 'multipole(P3) distribution';
        case 4
          directions = 'multipole(P4) distribution';
        case 7
          directions = 'multipole(P7) distribution';
        otherwise % NOT IMPLEMENTED
          directions = 'unknown distribution';
      end
    else
      % Guess that some 2D- or 3D- VMI-distribution with discrete polar angles (not nice for showing VMI image)
      directions = sprintf('%d different polar angles', size(f.rays,1));
    end
    distr_str = strrep(strrep(strrep(directions, 'multipole(',''), ') distribution',''), 'dipole(P1','dipole');
    directions = strrep(directions,' distribution',' distr');
    
    if ~ischar(source_point_pattern) && size(source_point_pattern,2) == 3 && size(source_point_pattern,1) == 1
      source_str = sprintf('FWHM=%.2f,%.2f,%.2f mm', 2*sqrt(2*log(2))*source_point_pattern(1), 2*sqrt(2*log(2))*source_point_pattern(2), 2*sqrt(2*log(2))*source_point_pattern(3));

    elseif length(source_point_pattern) == 1 && mod(source_point_pattern,1) == 0
      source_str = sprintf('Pattern %d', source_point_pattern); % TODO would need the value for source_point_spacing too
      if source_point_offset_vertical ~= 0
        error('source_point_offset_vertical is only supported when 3D Gaussian source point distribution by fly2-file.');
      end
    else
      source_str = mat2str(source_point_pattern);
    end
    if source_point_offset_y ~= 0 || source_point_offset_vertical ~= 0
      % (Offset along optical axis (y in SIMION & Matlab, horizontal in lab),
      %  Offset along vertical axis (z in SIMION, x in Matlab, vertical in lab) ,
      %  [if present: Offset along spectrometer (TOF) axis (x in SIMION, z in Matlab, horizontal in lab) )
      source_str = sprintf('%s, offset (%+.4g,%+.4g) mm', source_str, source_point_offset_y, source_point_offset_vertical);
    end
    if ~isempty(source_point_offset_z) && isfinite(source_point_offset_z)
      if source_point_offset_y ~= 0 || source_point_offset_vertical ~= 0
        % Append after (optical,vetical) offsets
        source_str = sprintf('%s,%+.4g) mm', strrep(source_str,') mm',''), source_point_offset_z);
      else
        % Just z-offset
        source_str = sprintf('%s, offset (0,0,%+.4g) mm', source_str, source_point_offset_z);
      end
    elseif get_argument_from_cells(variables, 'source_point_z')
      source_str = sprintf('%s,z=%.1fmm', source_str, get_argument_from_cells(variables, 'source_point_z', 'num'));
    end
    
    if length(energies) <= 7
      if all(energies < 0.1) %[eV]
        % All in meV
        energies_str = ['[', texformat_g10(energies/1E-3), '] meV'];
      else
        energies_str = ['[', texformat_g10(energies), '] eV'];
        energies_str = strrep(energies_str, '\cdot10^{-3}', 'm'); % use SI-prefix to shorten frequently occuring magnitude
      end
    else
      energies_str = ['[', strrep([texformat_g10(energies(1:3)) ',\_\_,' texformat_g10(energies(end-1:end))], ', ', ','), '] eV'];
    end
    if get_argument_from_cells(variables, 'initial_v_vertical', 'num', 0) ~= 0
      % Show vertical jet velocity with sign flipped to be positive for normal downwards direction
      energies_str = sprintf('%s + jet %.1fkm/s', energies_str, -get_argument_from_cells(variables, 'initial_v_vertical', 'num'));
    end
    title(sprintf('%s for %s, %s, %d particles/energy.\nIonization volume: %s. %s.', ...
      detector_name, energies_str, directions, min([max_shown_points, f.source.ion_count]), ...
      source_str, strrep(pot_str,'_','\_') ), 'FontSize',10);
    xlabel('Horizontal coordinate (optical axis) / mm')
    ylabel('Vertical coordinate (polarization axis) / mm')
    energies = fliplr(energies); % show original order for figure 4
    
    % Spectrum
    figure(5); clf;
    %set(5,'Position',[30 70 579 319]);
    set(5,'Position',[30 70 579 600]);
    subplot(3,1,1:2)
    if image_bin_width < 0.1 % using really small size, like camer pixels. Then keep pixel width of radial diagram
      radial_bin_width = image_bin_width;
    else % using a bit larger pixels so simulated 2D image looks OK even when not a huge amount of particles
      radial_bin_width = image_bin_width/2;
    end
    hist(radii(radii < halfwidth), (0:image_bin_max)*radial_bin_width);
    % xlim([20 30]); % zoom in for the 19&20 eV case
    xlim([0 35]); % fixed range
    %ylim([0 0.2*n]); % fixed range
    if any(original_energies < 7)
      ylim([0 1.56*radial_bin_width*n]); % fixed range based on bin width and particle count, not saturated by pure dipole in narrow bins
    else
      ylim([0 1.3*radial_bin_width*n]); % fixed range based on bin width and particle count, not saturated by pure dipole in narrow bins
    end
    xlabel('Radial coordinate / mm');
    ylabel('Electron count');
    
    % Energy calibration
    r = plane_radii_and_energies(:,1);
    K = plane_radii_and_energies(:,2); % NOTE: these are not discrete fixed energies but the projection to xy plane (thus a spread dependent on angular width)
    %which = K>0;
    which = K>0 & r>=1; % exclude very low radii, to avoid division by zero
    %which = K>0 & r>3; % exclude lowest radii if exlcuding zero energy would not be enough, to avoid division by zero
    r_uniform = 0:min([0.5, 5*radial_bin_width]):(max(r)+3*min([0.25, radial_bin_width]));
    kinetic_energy_calibration_2 = mean(K(which)./r(which).^2); % [eV/mm^2] averaging in energy (radius^2) domain, weights high radii more which we usually care more about
    kinetic_energy_calibration_1 = mean(sqrt(K(which))./r(which)).^2; % [eV/mm^2] averaging in velocity (radial) domain, reasonable weighting since since position measurement error independent of energy
    kinetic_energy_calibration_separate = NaN(size(energies));
    for ei = 1:length(energies)
      %kinetic_energy_calibration_separate(ei) = mean(K(which & K==energies(ei))./r(which & K==energies(ei)).^2); % [eV/mm^2] averaging in energy (radius^2) domain
      %kinetic_energy_calibration_separate(ei) = mean(sqrt(K(which & K==energies(ei)))./r(which & K==energies(ei))).^2; % [eV/mm^2] averaging in velocity (radial) domain
      this_energy = energies(ei);
      if length(energies) > 1
        separation_to_nearest = abs(energies - this_energy);
        separation_to_nearest = min(separation_to_nearest(separation_to_nearest > 0));
        which2 = abs(K-energies(ei)) < 0.5*separation_to_nearest; % find the trajectories for which this is the closest energy (actually the plane-projection can only make them lower, but allow for rounding errors and check symmetrically)
      else
        % When only one energy is simulated, not really a need to filter further. Of course the peak shape is asymmetric when not Abel inverted,
        % so a calibration using all trajectories will express some kind of adjusted caliration that makes the non-inverted peak centered on good enerty.
        % Perhaps better to take only the upper half of the radii, to at least come closer to the Abel-inverted ideal?
        which2 = K >= median(K);
      end
      kinetic_energy_calibration_separate(ei) = mean(sqrt(K(which2))./r(which2)).^2; % [eV/mm^2] averaging in velocity (radial) domain
    end
    kinetic_energy_calibration = [kinetic_energy_calibration_1 kinetic_energy_calibration_2 kinetic_energy_calibration_separate];
    kinetic_energy_calibration = sort(kinetic_energy_calibration);
    kinetic_energy_calibration(isnan(kinetic_energy_calibration)) = []; % skip NaN values (happens e.g. when just small radii)
    if length(kinetic_energy_calibration) > 2 % (unless NaN values were skipped already there should be three values)
      kinetic_energy_calibration([1 end]) = []; % skip lowest and highest value (the near-zero energy usually gives quite low-deviating calibration)
    end
    kinetic_energy_calibration = mean(kinetic_energy_calibration); %average the remaining
    
    source_str_optional = '';
    if length(energies) == 1
      source_str_optional = [', ' source_str];
    end
    title(sprintf('%dD %s %s, %s%s %s,\n%d particles each. %.3g mm bins. Vxy calibration = %.3fE-3 eV/mm^2 ', ...
      f.imaging_dimensions, detector_name, particle, directions, source_str_optional, energies_str, ...
      f.source.ion_count, radial_bin_width, kinetic_energy_calibration/1E-3), 'FontSize',10);
    % Not showing the used angular_halfwidth_selection_for_VMI, but its influence is not so big now that the initial velocity is projected onto plane.
    
    subplot(3,1,3); cla
%     plot(r.^2, K, '.', r_uniform.^2,r_uniform.^2*kinetic_energy_calibration,'--r');
%     xlabel('r^2 [mm^2]')
%     ylabel('Kinetic energy [eV]')
    plot(r(which).^2,K(which)./r(which).^2/0.001,'.', ...
         r_uniform.^2,ones(size(r_uniform))*kinetic_energy_calibration_2/0.001,':b', ...
         r_uniform.^2,ones(size(r_uniform))*kinetic_energy_calibration_1/0.001,'--b', ...
         energies./kinetic_energy_calibration_separate,kinetic_energy_calibration_separate/0.001,'or--', ...
         r_uniform.^2,ones(size(r_uniform))*kinetic_energy_calibration/0.001,'-r' );
    ylabel('Energy calibration [meV/mm^2]');
    if ~isnan(kinetic_energy_calibration)
      ylim([0.95 1.05] * kinetic_energy_calibration/0.001)
    end
    if max(kinetic_energy_calibration_separate)/0.001 > ylim()*[0;1]
      ylim([ylim()*[1;0], max(kinetic_energy_calibration_separate)/0.001])
    end
    if median(kinetic_energy_calibration_separate)/0.001 < ylim()*[1;0]
      ylim([median(kinetic_energy_calibration_separate)/0.001, ylim()*[1;1]])
    end
    if ~isempty(r_uniform)
      xlim([0 r_uniform(end).^2])
    end
    xlabel('r^2 [mm^2]')
    legend('Sim.', '2-norm', '1-norm', sprintf('1-norm tol.\nenergywise'), 'Combined', 'Location','EO')
    
    
    % Magnification etc. about mapping from initial spatial to detected spatial
    figure(7); clf
    set(gcf, 'Position', [379   123   560   821]);
    %radial_reference_size = 20; % [mm] Assuming ion detector ("B") is 40 mm diameter, 20 mm radius in all models.
    %radial_reference_size = 20/2; % [mm] Use half the detector radius as reference, to get ratio between third and first order coefficient of more intuitively magnitude
    radial_reference_size = 5; % [mm] Use a size within shown range, near achievable FWHM along optical axis
    if size(start_and_final_coords_where_nearly_axial_velocity,1) >= 5
      subplot(4,3,4:9);
      % Start with the lower, separate diagrams for the two axes, that clarify that for reflectron they can be quite different, 
      % inverted (negative first order) and with higher order terms (e.g. third order).
      halfrange = ceil(max(max(abs(start_and_final_coords_where_nearly_axial_velocity(:,1:2))))/0.5)*0.5; %[mm]
      % halfrange = 6; %[mm] fixed
      halfrange = 8; %[mm] fixed to include all that may reach detector with typical voltages
      x = -halfrange:0.1:halfrange;
      r_start = sqrt(sum(start_and_final_coords_where_nearly_axial_velocity(:,1:2).^2,2));
      r_det = sqrt(sum(start_and_final_coords_where_nearly_axial_velocity(:,3:4).^2,2));
      plot(r_start, r_det, '.g');
      hold on;
      %magnification_radial = r_start(r_start>1) \ r_det(r_start>1); %use points outside 1 mm radius, and do least-squares fit of linear coefficient (no constant term)
      magnification_radial = r_start \ r_det; % do least-squares fit of linear coefficient (no constant term)
      %pol_r = [magnification_radial 0];
      pol_r = polyfit(r_start, r_det, 1); % allow a non-centered output offset, without letting that affect magnification (slope) value
      %magnification_radial = pol_r(1);
      h_radial = plot([0 halfrange], magnification_radial*[0 halfrange], '-.g');
      plot([0 halfrange], pol_r(end) + pol_r(1)*[0 halfrange], '--g');
      % start_and_final_coords_where_nearly_axial_velocity [mm] contains the 4 columns: horizontal&vertical start, horizontal&vertical detected
      plot(start_and_final_coords_where_nearly_axial_velocity(:,1), start_and_final_coords_where_nearly_axial_velocity(:,3), '.k')
      % By fitting the polynomial after dividing position by radial_reference_size the coeffients become dimensionless,
      % and the relative weight between linear and third-order can be compared (corresponding to their relative influence at detector edge so even a third order of 1 is notable).
      pol_horizontal = polyfit(start_and_final_coords_where_nearly_axial_velocity(:,1)/radial_reference_size, start_and_final_coords_where_nearly_axial_velocity(:,3)/radial_reference_size, 3);
      h_fit_horizontal = plot(x, polyval(pol_horizontal, x/radial_reference_size)*radial_reference_size, '--k');
      plot(start_and_final_coords_where_nearly_axial_velocity(:,2), start_and_final_coords_where_nearly_axial_velocity(:,4), '.b')
      pol_vertical_order = 3;
      if std(start_and_final_coords_where_nearly_axial_velocity(:,2))*2*sqrt(2*log(2)) < 1.1
        % If the vertical source standard deviation is at most 1 mm (with some tolerance)
        % it is not meaningful to fit a third order coefficient, just random disturbance in the readout of dominant linear coefficient.
        pol_vertical_order = 2;
      end
      pol_vertical = polyfit(start_and_final_coords_where_nearly_axial_velocity(:,2)/radial_reference_size, start_and_final_coords_where_nearly_axial_velocity(:,4)/radial_reference_size, pol_vertical_order);
      h_fit_vertical = plot(x, polyval(pol_vertical, x/radial_reference_size)*radial_reference_size, '--b');
      if pol_vertical_order < 3
        pol_vertical = [NaN(1,3-pol_vertical_order) pol_vertical];
      end
      xlabel('Starting coordinate [mm]');
      ylabel('Imaged on detector [mm]');
      if image_from_detector_B
        axis([-halfrange halfrange -20 20]); % assume ion detector diameter 40 mm
      else
        axis([-halfrange halfrange -35 35]); % assume ion detector diameter 70 mm
      end
      horiz_vert_std = std(start_and_final_coords_where_nearly_axial_velocity(:,3:4)); %[mm]
      legend([h_fit_horizontal, h_fit_vertical, h_radial], ...
            sprintf('Horizontal centre %+.1f mm; linear %.1f, scaled quadratic %.1f & cubic %.1f.', pol_horizontal(end)*radial_reference_size, pol_horizontal(end-1), pol_horizontal(end-2), pol_horizontal(end-3)), ...
            sprintf('Vertical centre %+.1f mm; linear %.1f, scaled quadratic %.1f & cubic %.1f.', pol_vertical(end)*radial_reference_size, pol_vertical(end-1), pol_vertical(end-2), pol_vertical(end-3)), ...
            sprintf('Radial magnification %.2f. Alt. linear %+.2f, beyond %.1f mm offset', magnification_radial, pol_r(1), pol_r(end)), ...
            'Location','SO');
      title(sprintf('  Spatial mapping: Std. H=%.1f mm, V=%.1f mm. %.0f-mm-scaled cubic ratio H=%+.1f, V=%s', ...
              horiz_vert_std(1), horiz_vert_std(2), radial_reference_size, ...
              pol_horizontal(end-3)/pol_horizontal(end-1), ...
              iif(pol_vertical_order >= 3, sprintf('%+.1f', pol_vertical(end-3)/pol_vertical(end-1)), 'not fit') ...
              ), 'FontSize',9);
      set(gca, 'Position', [0.13 0.43 0.775 0.2419]); % move down a bit to not overlap other
    else
      pol_horizontal = NaN(1, 4);
      pol_vertical = NaN(1, 3);
      magnification_radial = NaN;
      mass_overlap_fraction = NaN;
      horiz_vert_std = NaN(1, 2);
      t_std = NaN;
      pol_t = NaN(1, 3);
    end
    
    if ~isempty(TOF_and_start_z0_all)
      % Show time of flight, and possibly a TOF histogram to see whether two test masses are resolved
      subplot(4,3,10:11);
      if isfield(f.source,'normal_source_point_z') && f.source.source_point_z ~= f.source.normal_source_point_z
        dz0_halfrange = max(abs(TOF_and_start_z0_all(:,2) - (f.source.source_point_z - f.source.normal_source_point_z)));
        dz0 = -dz0_halfrange:0.01:dz0_halfrange;
        % Shift the shown range to be centered on the simulated source points,
        % while keeping 0 of the axis to mean the normal_source_point_z.
        dz0 = dz0 + (f.source.source_point_z - f.source.normal_source_point_z);
      else
        dz0_halfrange = max(abs(TOF_and_start_z0_all(:,2)));
        dz0 = -dz0_halfrange:0.01:dz0_halfrange;
      end
      is_heaviest_mass = TOF_and_start_z0_all(:,3)==max(TOF_and_start_z0_all(:,3));
      plot(TOF_and_start_z0_all(:,2), TOF_and_start_z0_all(:,1)*(1E-9/1E-6), '.m');
      xlabel('Starting coordinate along spectrometer axis [mm]');
      axis tight
      yl = quantile(TOF_and_start_z0_all(:,1), [0.02 0.98]) * (1E-9/1E-6) + [-0.01 0.01];
      yl = mean([ylim(); yl]);
      hold on
      pol_t = polyfit(TOF_and_start_z0_all(is_heaviest_mass,2), TOF_and_start_z0_all(is_heaviest_mass,1), 3);
      plot(dz0, polyval(pol_t, dz0)*(1E-9/1E-6), '--r');
      xlim(dz0([1 end]));
      ylim(yl);
      t_std = std(TOF_and_start_z0_all(is_heaviest_mass,1));
      ylabel(['TOF [' 956 's] for ' sprintf('%s u/e', strrep(mat2str(unique(TOF_and_start_z0_all(:,3))'), ' ', ','))]); % 956 is the Unicode for ? (upright mu). It seems a plaintext ? sometimes get lost when saving/loading/running the file (not by F9)
      title(sprintf('{\\itt} / ns = %.0f %+.1f ({\\itz}_0/mm) %+.1f ({\\itz}_0/mm)^2 %+.1f ({\\itz}_0/mm)^3      ', ...
              pol_t(end), pol_t(end-1), pol_t(end-2), pol_t(end-3)), 'FontSize',8.5);
      set(gca, 'YGrid','on', 'TickLength',[0.014 0.02]);
      set(gca, 'Position', [0.13 0.08 0.55 0.1577]); % move down a bit to not overlap other, and widen
      subplot(4,3,12);
      % Small histogram to represent the TOF-peak
      t_hist_bins = floor(yl(1)*1E3):2:ceil(yl(2)*1E3); %[ns] use 2 ns bins like the 500 MHz sampling of waveform ADC card
      [t_hist, t_hist_bins] = hist(TOF_and_start_z0_all(:,1), t_hist_bins);
      barh(t_hist_bins*(1E-9/1E-6),t_hist, 1, 'EdgeColor','none');
      ylim(yl);
      xlim([0 max(t_hist)+1]);
      set(gca, 'YGrid','on', 'YTickLabel',[], 'TickLength',[0.025 0.025]);
      set(gca, 'Position', [0.6916 0.08 0.2134 0.1577]); % move down a bit to not overlap other
      TOF_count = sum(is_heaviest_mass);
      if abs(sum(~is_heaviest_mass)-TOF_count) < 0.05*TOF_count % if there is a lower test mass (with a similar number of non-miss trajectories)
        TOF_count = mean([TOF_count sum(~is_heaviest_mass)]); % use the average count to define the fraction
        % Consider the start time of the heavier mass to be its 0.5 percent percentile, minus half of the 2 ns sampling time of our acquisition card
        t_start_of_heavy = quantile(TOF_and_start_z0_all(is_heaviest_mass,1), 0.005) - 1; %[ns]
        mass_overlap_fraction = sum(TOF_and_start_z0_all(~is_heaviest_mass,1) > t_start_of_heavy) / TOF_count;
      else
        mass_overlap_fraction = NaN;
      end
      title(sprintf('    Std. = %.1f ns = {\\itt}/%s\n    %.0f%% overlap @%.0f u', ...
              t_std, sprintf('%.0f', mean(TOF_and_start_z0_all(:,1))/t_std), ...
              100*mass_overlap_fraction, max(TOF_and_start_z0_all(is_heaviest_mass,3))), ... 
            'FontSize',9);
      xlabel('Count');
    end    
    
    subplot(4,3,1:3); % Magnification
    magnification_avg = mean(magnifications_and_start_radii(magnifications_and_start_radii(:,2) > 2,1)); % average for starting points further than 2 mm from centre (to not divide by near zero)
    magnification_RMS = rms(magnifications_and_start_radii(magnifications_and_start_radii(:,2) > 2,1)); % RMS (so that also negative magnification contributes to VMI-focusing imperfection)
    plot(magnifications_and_start_radii(:,2), magnifications_and_start_radii(:,1), '.k', ...
      [0 10],[1 1]*magnification_avg,'--b', [0 10],[1 1]*magnification_RMS,'-r' );
    xlabel('Starting radial coordinate [mm]')
    ylabel(sprintf('Spatial-mode magnification\n({\\bfimaged}\\cdot{\\bfsource})/|{\\bfsource}|^2'), 'FontSize',10)
    %legend('Sim.', 'Avg.', 'RMS.')
    if magnification_RMS > 0.5 % seems like spatial mode
      if image_from_detector_B
        % Use wide range for the reflectron detector, where magnification is often negative too (due to mirror or focusing?)
        ylim([-20 20])
      else
        % Normal spatial mode on VMI-detector. Use fixed axis range so slope can be judged
        ylim([0 2])
      end
    end
    title(sprintf('%dD %s %s, %s %s,\nFor max %.0f\\circ: spatial magnification avg.=%s, RMS=%.3G', ...
      f.imaging_dimensions, detector_name, particle, directions, energies_str, ...
      angular_halfwidth_selection_for_magnification, texformat_g10(magnification_avg,3), magnification_RMS), ...
      'FontSize',10);
    % TODO IMPROVEMENT: could subtract offset due to velocity, using energy calib. for more precise magnification, even without tight angular filter
    
    % Append to log:
                            % workbench_name, job_name, offsets, potential_offset_str, potential_offset_name, potential_offset_amount, potentials_orig, potentials_offset, potentials, variables
    reflectron_log_inputs = [reflectron_log_inputs; {workbench_name, job_name, offset_list(offset_index,:), potential_offset_str, potential_offset_name, potential_offset_amount, potentials_orig, potential_offset, potentials, variables}];
    reflectron_log_results = [reflectron_log_results; ...
        ... % mean(H,V,t), std(H,V,t), [H.linear, H.cubic, V.linear, magnification_radial], [t/std(t), mass overlap[%], dT/dz0(linear term)]
        [pol_horizontal(end)*radial_reference_size, pol_vertical(end)*radial_reference_size, mean(TOF_and_start_z0_all(:,1))], ... % [mean(H,V,t)], 
        [horiz_vert_std(1), horiz_vert_std(2), t_std], ... % [std(H,V,t)], 
        [pol_horizontal(end-1), pol_horizontal(end-3), pol_vertical(end-1), magnification_radial], ... % H.linear, H.cubic, V.linear, magnification_radial
        [mean(TOF_and_start_z0_all(:,1))/t_std, 100*mass_overlap_fraction,  pol_t(end-1)] % t/std(t), mass overlap[%], dT/dz0(linear term)
      ];
    
    
    % Show more like camera image
    figure(4); clf;
    set(4, 'Position', [900 59 798 698]);
    imagesc8px([-halfwidth halfwidth], [-halfwidth halfwidth], image_histogram, image_bin_width);
    colorbar
    daspect([1 1 1]);
    if not(show_starts)
      grid on
    end
    set(gca,'GridColor', [1 1 1]*0, 'GridAlpha',0.5);
    % Colouring
    clim = get(gca,'CLim');
    if ~show_starts
      if any(energies < 0.03)
        redundancy_factor = sum(energies <= min(energies)+1E-4);
        redundancy_factor = inf; % OPTION to force auto-scaled intensity range (cropping at a fixed fraction of max)
        if redundancy_factor > 4 % (actually the redundancy factor seems to work in currently used particle lists, so this case is not important)
          % Multiple "particles" (i.e. mass and energy combinations) probably overlap,
          % particularly for the case of low energy ions this means the scaling by source.ion_count
          % is not enough to guess a good intensity range.
          % Don't try to use fixed intensity range, let it autoscale, and then crop to 2/3 of max
          set(gca,'CLim', round(clim .* 2/3));
        else
          % Use a fixed intensity range, estimated heuristically based on the 
          % number of particles and image histogram bin width.
      %   set(gca,'CLim', floor(clim .* 0.005));
      %  set(gca,'CLim', floor(clim .* (0.005 / image_bin_width)));
      %  set(gca,'CLim', [0 flys(particle_indices(1)).source.ion_count * image_bin_width.^2 * 0.05]);
      %  set(gca,'CLim', [0 flys(particle_indices(1)).source.ion_count * image_bin_width.^2 * 0.1]);
          set(gca,'CLim', [0 flys(particle_indices(1)).source.ion_count * image_bin_width.^2 * 0.08 * redundancy_factor]);
        end
        
      %   set(gca,'CLim',[0 100]); % when 0eV is in simulation
      elseif length(particle_indices) == 5
         % when innermost ring has much lower kinetic energy, 
        set(gca,'CLim', round(clim .* 0.55));
      end
      
    else % show_starts: 
      % Draw a ring for the extractor electrode's inner diameter
      h = rectangle(gca, 'Position',[-10 -10 20 20 ], 'Curvature',[1 1], 'EdgeColor','k', 'LineStyle',':');
    end
    
    if exist('custom_colormap_jet_like','file')
      colormap(custom_colormap_jet_like())
    else
      % colormap(jet(128)); brighten(-0.3);
      colormap(jet(128));
      % colormap(flipud(gray(128))); brighten(-0.2);
      % colormap(flipud(gray(128))); brighten(-0.35);
    end
    %title(sprintf('%s for %s, %s, %d particles/energy. Pixels of %.3g mm.\nIonization %s. %s.', detector_name, energies_str, directions, f.source.ion_count, image_bin_width, source_str, strrep(pot_str,'_','\_')), 'Interpreter','TeX');
    title(sprintf('%s for %s, %s, %s, %d/energy, %.1f %% missed. Pixels of %.3g mm.\nIonization %s. %s.', ...
      detector_name, energies_str, particle, distr_str, f.source.ion_count, 100*miss_fraction, ...
      image_bin_width, source_str, strrep(pot_str,'_','\_')), 'Interpreter','TeX', 'FontSize',10);
    if show_starts
      title(sprintf('Starting points, %d particles/energy. Pixels of %.3g mm.\nIonization %s', ...
        f.source.ion_count, image_bin_width, strrep(pot_str,'_','\_')));
    end
    if image_from_detector_B
      set(gca,'XDir','reverse')
      % coords(:,1), increases towards the optical exit side (assuming small angle of detector B), 
      % i.e. the usual coordinate system as seen from VMI side (detector A).
      % However, the optical mirror inside chamber used for seeing detector B adds a
      % horizontal mirroring, and to ease the comparison with experimental images,
      % we reproduce that mirroring here. Thus right in the simulated image will be right.
      % Keeping the coordinates called negative on right side, so the sign of magnification coefficients etc. is not affected.
      xlabel('{\itNegative on right to appear as through mirror:} Horizontal coordinate (optical axis) / mm')
    else
      xlabel('Horizontal coordinate (optical axis) / mm')
    end
    ylabel('Vertical coordinate (polarization axis) / mm')
    
    filename = sprintf('FWHM%.1f,%.1fmm@%.4g,%.1fmm', ...
          2*sqrt(2*log(2))*source_point_pattern(2), 2*sqrt(2*log(2))*source_point_pattern(3), ...
          source_point_offset_y, source_point_offset_vertical);
    if isfield(f.source,'normal_source_point_z')
      % Verify that correct offset
      if f.source.source_point_z ~= f.source.normal_source_point_z
        if ~isempty(source_point_offset_z) && source_point_offset_z ~= f.source.source_point_z - f.source.normal_source_point_z
          warning('Unesual source_point_offset_z. Parameter %.2f mm but used %.2f mm', source_point_offset_z, f.source.source_point_z - f.source.normal_source_point_z);
        end
        source_point_offset_z = f.source.source_point_z - f.source.normal_source_point_z;
      end
    end
    if ~isempty(source_point_offset_z) && isfinite(source_point_offset_z)
      if strcmp(filename(end-1:end), 'mm')
        filename = filename(1:end-2);
      end
      filename = sprintf('%s,%.1fmm', filename, source_point_offset_z);
    elseif get_argument_from_cells(variables, 'source_point_z')
      filename = sprintf('%s,z%.1fmm', filename, get_argument_from_cells(variables, 'source_point_z', 'num'));
    end
    % Optional shortening and tricks
    filename = strrep(filename, 'FWHM', 'W');
    filename = strrep(filename, 'mm@', '@');
    filename = strrep(filename, '-', ' -'); % by putting a space before minus signs, negative values will get sorted before zero and positive in Windows Explorer (order among the negative if multiple negative has not been checked)
    if get_argument_from_cells(variables, 'initial_v_vertical')
      initial_v_vertical = get_argument_from_cells(variables, 'initial_v_vertical', 'num');
      if initial_v_vertical == 0
        filename = sprintf('%s,jet%.1f', filename, initial_v_vertical); % 0 without + or - sign
      elseif initial_v_vertical < 0
        filename = sprintf('%s,jet%.1f', filename, -initial_v_vertical); %[km/s] without the minus sign for normal downwards direction
      else
        filename = sprintf('%s,jet%+.1f', filename, initial_v_vertical); %[km/s] with a plus sign for unusual upwards direction
      end
    end
    if show_starts
      filename = sprintf('source %s', filename);
    else % show voltage info in filenames
      filename = sprintf('%s %s, %s', distr_str, pot_str_w_offset, filename);
    end
    if all(mass>=1) && length(mass) == 1 % OPTION to automatically enable mass in filename
      filename = sprintf('%s %su', filename, mat2str(unique(mass))); % OPTION include mass
      filename = strrep(filename, '000u', 'ku'); % shortne to kDa mass unit (shown as 'ku' still)
    end
    
    %if true | size(names_workbenches_and_potentials,1) > 2 % OPTION to include an index to distinguish similar potential lists
    %  filename = sprintf('j%03d %s', job+first_job-1, strrep(strrep(filename, [distr_str ' '], ''), [workbench_name ' '], ''));
    %end
    if multipole ~= 0
      % OPTION: skip distr_str and workbench_name to shorten filenames when including lots of other info in filename
      filename = strrep(strrep(filename, [distr_str ' '], ''), [workbench_name ' '], '');
      filename = strrep(filename, 'distribution', 'distr');
    else
      filename = strrep(filename, 'isotropic distribution', 'isotr');
    end
    filename = [job_name ' ' filename]; % OPTION include job_name
    
    %export_fig_multiformat(3,['dots ' filename], 'PNG')%,'PDF')
    export_fig_multiformat(4,['image ' filename], 'PNG')%, 'FIG')%,'PDF')
    if ~image_from_detector_B
      export_fig_multiformat(5,['spectrum ' filename], 'PNG')
    end
    export_fig_multiformat(7,['magn ' filename], 'PNG') % magnification (spatial mapping, and also z0-to-TOF and mass resolution)
    % OPTION to save also the numeric data:
%     save(['sim ' strrep(filename,'.0mm','mm') '.mat'], 'variables', 'workbench_name', 'mass_in_u', 'original_energies', 'energies', 'source_point_pattern', 'source_point_offset_y','source_point_offset_vertical','detector_position','potentials','pot_str','multipole','n','directions','show_starts','image_bin_width','flys');
    % save(['histogram ' strrep(filename,'.0mm','mm') '.mat'], 'variables', 'workbench_name', 'mass_in_u', 'original_energies', 'energies', 'source_point_pattern', 'source_point_offset_y','source_point_offset_vertical','detector_position','potentials','pot_str','multipole','n','directions','show_starts','image_bin_width','image_histogram');


    %% End of one iteration
  end % end of dV_index-loop
  end % end of dV_seq_index-loop
  % Increment the combinatorial_indices, like a positional number in base length(potential_offset_values) (but each 1-based for Matlab array indexing rather than 0-based)
  for c = length(combinatorial_indices):-1:1
    % Increment at this position.
    combinatorial_indices(c) = combinatorial_indices(c) + 1;
    if combinatorial_indices(c) <= length(potential_offset_values)
      % Still a valid index at this position. Don't need to increment preceeding positions now.
      break
    else
      % Last index (for potential_offset_values) has been used at this position.
      % Reset to 1 and proceed to increment the preceeding position.
      combinatorial_indices(c) = 1;
    end
  end
  end % end of dV_combinatorial_index 
  
  end % end offset_index-loop
end % end job-loop
if ~isempty(potential_offset_components_sequential) || ~isempty(potential_offset_components_combinatorial)
  save('log tmp.mat', 'reflectron_log_inputs', 'potential_offset_components', 'potential_offset_values', 'reflectron_log_results');
end

if ~isempty(potential_offset_components_sequential) && ~isempty(potential_offset_values)
  % Show effect of each component
  if length(potential_offset_values) == 2
    % Take the differential
    reflectron_log_diffs = reflectron_log_results(2:2:end,:) - reflectron_log_results(1:2:end,:);
    reflectron_log_magnitude = max(abs(reflectron_log_results(2:2:end,:)), abs(reflectron_log_results(1:2:end,:))); % find something nonzero to divide with
    reflectron_log_diff_ratios = reflectron_log_diffs ./ reflectron_log_magnitude;
    N = size(reflectron_log_diffs,1) / size(potential_offset_components_sequential,1); % the number of jobs (counting position offset as separate jobs)
    grouped_diffs      = NaN(N, size(potential_offset_components_sequential,1), size(reflectron_log_results,2));
    grouped_magnitudes = NaN(N, size(potential_offset_components_sequential,1), size(reflectron_log_results,2));
    grouped_ratios     = NaN(N, size(potential_offset_components_sequential,1), size(reflectron_log_results,2));
    for c = 1:size(potential_offset_components_sequential,1)
      grouped_diffs(:,c,:)      = reflectron_log_diffs(c:size(potential_offset_components_sequential,1):end,:);
      grouped_magnitudes(:,c,:) = reflectron_log_magnitude(c:size(potential_offset_components_sequential,1):end,:);
      grouped_ratios(:,c,:)     = reflectron_log_diff_ratios(c:size(potential_offset_components_sequential,1):end,:);
    end
%     grouped_diffs(:,:,3)       = grouped_diffs(:,:,3) / 1E3; % us rather than ns for absolute TOF
%     grouped_magnitudes(:,:,11) = grouped_magnitudes(:,11) / 1E3; % kDa rather than u for t/std(t))
    result_magnitude = transpose(squeeze(mean(mean(grouped_magnitudes,2),1))); % average for component and job(voltage setting), to get one magnitude per result column
    % The ratios are relative differences: change/typical_value (to compare all results on same dimensionless scale), not divided by the voltage difference that caused it
    ratios = grouped_diffs ./ repmat(reshape(result_magnitude,1,1,length(result_magnitude)), size(grouped_diffs,1), size(grouped_diffs,2), 1);
    avg_ratios = squeeze(mean(ratios,1)); % average the jobs (voltage settings), indexed by (potential offset component, result column)
    variation = squeeze(std(ratios,[],1)) ./ squeeze(mean(abs(ratios),1));  % (std/mean(abs)), indexed by (potential offset component, result column)
    consistent = variation < 1.2; % == std < 1.2*mean(abs), indexed by (potential offset component, result column)
    show = consistent & avg_ratios >= 5E-3; % the difference between test voltage restuls is more than 0.5% of typical magnitude
    %Columns of reflectron_log_results:  mean(H,V,t), std(H,V,t), [H.linear, H.cubic, V.linear, magnification_radial], [t/std(t), mass overlap[%], dT/dz0(linear term)]
    column_names = {'H' 'V' 'T','stdH' 'stdV' 'stdT', 'H.1' 'H.3' 'V.1' 'R.mg', 'm.T' 'm.o' 'WMcL'};
    for cm = 1:size(potential_offset_components_sequential,1)
      if any(show(cm,:))
        columns = find(show(cm,:));
        [~,reordering] = sort(avg_ratios(cm,columns));
        columns = columns(fliplr(reordering)); % sort in descending order of ratio
        if length(columns) > 2
          % If many columns, skip columns with a ratio below 7.5% and below 1/5 of max among columns
          columns(avg_ratios(cm,columns) < min([0.075 0.2*max(avg_ratios(cm,columns))])) = [];
        end
        fprintf('''%s''-component affects: ', potential_offset_names_sequential{cm});
        for col = columns
          fprintf('%4s %4.3g%%, ', column_names{col}, avg_ratios(cm,col)*100);
        end
        fprintf('\n');
      end
    end
    for col = 1:size(reflectron_log_results,2)
      if any(show(:,col))
        components = find(show(:,col));
        [~,reordering] = sort(avg_ratios(components,col));
        components = components(flipud(reordering)); % sort in descending order of ratio
        if length(components) > 2
          % If many columns, skip columns with a ratio below 5% and below 1/5 of max among columns
          components(avg_ratios(components,col) < min([0.05 0.2*max(avg_ratios(components,col))])) = [];
        end
        fprintf('%4s-result affected by: ', column_names{col});
        for cm = transpose(components)
          fprintf('''%s'' %4.3g%%, ', potential_offset_names_sequential{cm}, avg_ratios(cm,col)*100);
        end
        fprintf('\n');
      end
    end
  end
end