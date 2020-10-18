% Function to run voltage optimizations of one- or two-sided electrostatic spectrometers.
% 
% The launching and readout of results into a scalar or 2-dimensional score is handled
% by optimize_potentials_objective, and the optimization can be done with various algorithms
% including standard deterministic algorithms (e.g. 'active-set' or 'interior-point' 
% while 'trust-region-reflective' is too slow) or randomizing genetic-style optimization
% either 'genetic' returning a single result or 'pareto' returning a few results
% that represent different (all good in some sense) compromises between two scores
% (i.e. electron-VMI-score vs. ion-mass-resolution score in a two-sided spectrometer).
% 
% Instead of running for fixed range of potentials, start with manually (or rough automatic)
% selected setting and use nonlinear optimizer to improve it.
% PARAMETERS
%  workbench_name
%  potentials_initial  I-by-P (Default: [])
%                  If potentials_initial is empty, it will be set to the arithmentic mean of potentials_min & potentials_max.
%                  Normally a single row (I=1) is given, but for genetic algorithms additional rows will be used to create the initial population.
%  potentials_min
%  potentials_max
%                  The potentials_... specify an optional initial guess and ranges for the optimization.
%                  If min and max are equal for an index, that potential will be held fixed (subject to overriding by potential_update_matrix).
%                  If a genetic algorithm is used, potentials_initial is ignored.
%  potential_update_matrix  (Default: eye(length(potentials_min)))
%                  Allows some of the non-optimized ("fixed") potentials to instead be given
%                  by linear combination of other potentials.
%                  Example:
%                    potential_update_matrix = eye(length(potentials_min));
%                    potential_update_matrix(B_LAST,:) = (1:length(potentials_min)) == B_FIRST;
%                    potential_update_matrix(5,:) = [0 0.5 0 0 0 0.5];
%                  would make the electrode at index B_LAST use the potential of that with index B_FIRST,
%                  and would let the 5th electrode be the average of the 2nd and 6th electrode.
%
%  absolute_index (Default: [])
%                  Empty array to directly use the potentials as parameters
%                  An index (>=1, <=length(potentials_min)) to use a probably more efficient parametrization
%                  wher one parameter is the A-repeller potential (coefficient wrt. initial guess)
%                  and other parameters are corefficients wrt. A-repeller. Thus rescaling lens for
%                  higher/lower kinetic energy requires only one parameter to be modified.
%  quantity       (Default: 'y_rel mix')
%  algorithm      (Default: 'active-set', or 'genetic' if potentials_initial is empty)
%                 The main option are 'active-set', 'interior-point' or 'genetic'.
%                 Suffices of ' quick', ' fine', ' quick-debug' or ' wide' (and maybe other) can be used 
%                 to adjust options.
%  batch_name     (Default: '' to not write to log file)
%                 Name of log file to append result to.
%                 If batch_name is 'just evaluate initial simulation', no optimization will be run beyond the potentials_initial.
%
% RETURNS
%  potentials     Result of optimization.
%  lowest_cost    Result of optimization. Normally scalar but
%                 costs_before_averaging (which may be a 1-by-N array) is returned if batch_name is 'just evaluate initial simulation'.
%  output         Info from the optimizer: stop condition, iteration count etc.
%                 If batch_name is 'just evaluate initial simulation', this returns the flys struct instead.
%  str            Message (optionally printed to log) about the optimization parameters and result.
%                 If batch_name is 'just evaluate initial simulation', this is omitted.
%  dependent      Tells which potentials are linear combinations of other potentials.
%                 If batch_name is 'just evaluate initial simulation', this is omitted.
% 
% SEE ALSO
% optimize_potentials_objective.m, optimize_potentials_loop.m
% 
function [potentials, lowest_cost, output, str, dependent] = optimize_potentials(workbench_name, ...
  potentials_initial, potentials_min, potentials_max, potential_update_matrix, absolute_index, ...
  quantity, algorithm, batch_name)

if contains(workbench_name, 'both8')
  threads = 1 % Seems better when each SIMION instance requires too much memory (so swapping starts and loading takes several seconds)
else
  threads = 2 % number of threads (max 2 allowed on laptop 2016, with default profile...?)
  % threads = 1 % DEBUG/OPTION Seems better when each SIMION instance requires too much memory (so swapping starts and loading takes several seconds)
end

% global cached_recent
% cached_recent = []; % clear cache
% clar cached_recent % close access to the globlal variable


%% Further settings
quasiisotropic_for_VMI = ~isempty(strfind(quantity, 'iso')); % 2015-12-30 to use _isotropic==3

if contains(quantity, 'FWHM5') % added 2020-02-25
  % A bit larger 3D-Gaussian than default, still not very long along optical axis
  % Reasonable test case for both8eu, and some vertical width (it does affect reflectron transmittance/loss) to perhaps emulate jet temperature spread or the vertical dispersion by mass&TOF (due to jet velocity)
  source_point_distribution = [0.2 5 1] / (2*sqrt(2*log(2))); % convert from FWHM [mm] along (TOF axis, optical axis, vertical (jet) axis)
elseif contains(quantity, 'FWHM3') % added 2020-02-25
  % Intermediate between 5mm and the small default (unrealistically small but good to optimize on, often giveing OK result also with 5 mm FHWM)
  %previously: source_point_distribution = [0.2 3 1] / (2*sqrt(2*log(2))); % convert from FWHM [mm] along (TOF axis, optical axis, vertical (jet) axis)
  % 2020-03-06 since FWHM3 with 5kV still didn't give mass resolving powers above 1, maybe the other axes than y really do matter (since not Wiley-McLauren ideal)
  source_point_distribution = [0.12 3 0.5] / (2*sqrt(2*log(2))); % convert from FWHM [mm] along (TOF axis, optical axis, vertical (jet) axis)
elseif contains(quantity, 'FWHM8') % added 2020-02-25
  % Quite long along optical axis (a notable amount of ions will be blocked by first electrode's aperture)
  source_point_distribution = [0.2 8 0.5] / (2*sqrt(2*log(2))); % convert from FWHM [mm] along (TOF axis, optical axis, vertical (jet) axis)
else % as used 2016
  source_point_distribution = []; % use default 3D-Gaussian [0.06 1 0.06] in sigma, corresponding to [0.1413 2.3548 0.1413 mm FWHM
end


if ~isempty(strfind(quantity, 'mass'))
  % Mass resolution. Not intending to momantum-image the ions, just get width of TOF-distribution.

  % TODO: multiple masses
  %energies = [1E-4, 0.01 1 2, 3]; % [eV] -- see script_WMcL_attempts.m. To have >= 1eV is quite challenging since any miss is given high penalty.
  energies = [1E-4, 0.01 1]; % [eV] -- easier to start with, but perhaps not as relevant for dissociation which is the challenging case
  energies = [1]; % [eV] -- to avoid averaging  
%   mass = 100; % will get positive charge

  % To allow using real spacing to judge mass resolution, simulate neighbouring masses of same kinetic energy.
  % Since 1eV (aim 200u) was found most challenging  in 1D-analytical optimizations (among <meV, 0.01eV@300u,1eV@200u,2eV@100u,3eV@1u)
  % and somewhat more on M_dissociation(low mass, high energy) than M_recoil(high mass, low energy) side, it should be a good single test energy.
  % One could think 150u@1.5eV would be a relevant midpoint, but since the 1D-analytic work was using the perhaps optimistic dz=60um
  % it is reasonable to now optimize for a more |Q_z|-sensitive case (low energy) than to experimentally find that |Q_z| became too large.
  % energies = [1 1]; mass = [100 101]; % Used until 2015-04-28
%   energies = [1 1]; mass = [150 151]; % ?
  energies = [1 1]; mass = [199 200]; % ?Used from 2015-04-28, since 100u was achieved and 200u is the aim for 1eV.
  % Note: Optimising for 1eV may mean low angular collection for 3eV H+. TODO: test.
  
  particle_count = 500;
  
  % Allow suffixes on the form ' 0.001eV@200u' to set which energies and mass to use
  match = regexp(quantity, ' ([0-9.]+)eV@([0-9]+)u', 'tokens');
  if ~isempty(match)
    energies = []; mass = [];
    for i = 1:length(match) % m
      e = str2double(match{i}{1});
      m = str2double(match{i}{2});
      if ~isempty(e) && e >= 0 && e <= 100 && ~isempty(m) && m >= 5 && isfinite(m)
        energies = [energies, e e  ]; %[eV]
        mass     = [mass,     m m+1]; %[u]
      end
    end
    if length(mass)/2 >= 4
      particle_count = 300; % reduce the number of particles per mass somewhat, to speed up. An average of results for all masses will be used to drive optimization.
    elseif length(mass)/2 >= 2 % i.e. at least two groups of m,m+1
      particle_count = 300; % reduce the number of particles per mass somewhat, to speed up. An average of results for all masses will be used to drive optimization.
    end
    if isempty(mass)
      error('The quantity to optimize contained suffixes that could not be interpreted as a list of <energy>eV@<mass>u.')
    end
  end
  
  % Customization of which electron energies
  if ~isempty(strfind(quantity, ' 5en ')) 
    % Include also intermediate energies? I.e. 5 energies instead of 3?
    im = 1; lo = 1;
  elseif ~isempty(strfind(quantity, ' h3en')) 
    % Only the higher three energies among the 5 defined.
    im = 1; lo = [];
  else
    im = []; lo = 1;
  end
  
  if ~isempty(strfind(quantity, 'y_')) % append also low-energy electron VMI simulation
    if length(mass) == 1
      mass = repmat(mass,1,length(energies));
    end
    mass     = [mass    , 0*lo    0*im*lo     0    0*im     0]; % electron, will get negative charge
    energies = [energies, 0.4*lo, 1.95*im*lo, 3.5, 5.05*im, 6.6]; % [eV] 1.55*2 eV spacing. Skip intermediate energies to speed up this combined ion + electron optimization.
%     variables = {'_number_of_directions',8, 'particles_per_direction',30, 'trajectories',true}; % was too noisy with 20/dir.
    number_of_directions = 8; particle_count = 38; % was too noisy with 20/dir. To reach 38*8=304 (near ordinary 300) particles in total.
  elseif ~isempty(strfind(quantity, 'y10_')) % append also low-energy electron VMI simulation, up to 10 eV (added 2020-02-25)
    if length(mass) == 1
      mass = repmat(mass,1,length(energies));
    end
    mass     = [mass    , 0*lo    0*im*lo     0    0*im     0]; % electron, will get negative charge
    energies = [energies, 0.67*lo, 3.0*im*lo, 5.33, 7.67*im, 10.0]; % [eV] % rescaled from 30eV-list
%     variables = {'_number_of_directions',8, 'particles_per_direction',30, 'trajectories',true}; % was too noisy with 20/dir.
    number_of_directions = 8; particle_count = 38; % was too noisy with 20/dir. To reach 38*8=304 (near ordinary 300) particles in total.
  elseif ~isempty(strfind(quantity, 'y30_')) % append also high-energy electron VMI simulation
    if length(mass) == 1
      mass = repmat(mass,1,length(energies));
    end
    mass     = [mass    , 0*lo    0*im*lo     0    0*im     0]; % electron, will get negative charge
    energies = [energies, 2*lo, 9*im*lo, 16, 23*im, 30]; % [eV] 1.55*2 eV spacing. Skip intermediate energies to speed up this combined ion + electron optimization.
    %variables = {'_number_of_directions',8, 'particles_per_direction',30, 'trajectories',true}; % was too noisy with 20/dir.
    number_of_directions = 8; particle_count = 38; % was too noisy with 20/dir. To reach 38*8=304 (near ordinary 300) particles in total.
  elseif ~isempty(strfind(quantity, 'y20_')) % scale the y30 to intermediate energy range, 10-20 eV probably more relevant than 30eV.
    if length(mass) == 1
      mass = repmat(mass,1,length(energies));
    end
    mass     = [mass    , 0*lo    0*im*lo     0    0*im     0]; % electron, will get negative charge
    %energies = [energies, 1,  4.5*im, 8,  11.5*im, 15  ]; % [eV] 3.5*2 eV spacing. Skip intermediate energies to speed up this combined ion + electron optimization.
    energies = [energies, 1.5*lo, 6*im*lo, 10.5, 15*im,  19.5]; % [eV] 4.5*2 eV spacing. Skip intermediate energies to speed up this combined ion + electron optimization.
    number_of_directions = 8; particle_count = 38; % was too noisy with 20/dir. To reach 38*8=304 (near ordinary 300) particles in total.
  elseif ~isempty(strfind(quantity, 'y42_')) % append also high-energy electron VMI simulation
    if length(mass) == 1
      mass = repmat(mass,1,length(energies));
    end
    mass     = [mass    , 0*lo    0*im*lo     0    0*im     0]; % electron, will get negative charge
    energies = [energies, [2*lo, 9*im*lo, 16, 23*im, 30] *42/30]; % [eV] scale to higher energy
    number_of_directions = 8; particle_count = 38; % was too noisy with 20/dir. To reach 38*8=304 (near ordinary 300) particles in total.

  elseif quasiisotropic_for_VMI
    % Run only ions but using discrete directions to check VMI (not a main usage case)
    number_of_directions = 8; particle_count = 38; 
  else
    % When not including eVMI, and not deterministic approximately 3D-isotropic desired, retain old random 3D-isotropic
    number_of_directions = 1;
    %(number_of_directions = 2; would now give dipole distribution (along y))
  end
  
	if ~isempty(strfind(quantity, 'many'))
    % optionally ten times as many particles
    particle_count = 10 * particle_count;
  end
  variables = {'_number_of_directions',number_of_directions, 'particles_per_direction',particle_count, 'trajectories',true}; % 2015-12-30 when _isotropic==3 will be used always (regardless of whether y30_ is included or not)
  
  if quasiisotropic_for_VMI
    variables{end+1}='_isotropic'; variables{end+1}=3; % From 2015-12-30 new 3D-isotropic VMI-style
  elseif number_of_directions == 1
    variables{end+1}='_isotropic'; variables{end+1}=1; % Only mass resolution (no (e)VMI), as before 2015-12-30
  else
    % Use default distribution. For current workbenches (Lua program) this would be VMI with _isotropic=2.
  end

elseif ~isempty(strfind(quantity, 'spheres'))
  % Not for optimization, but to show a fine simulated VMI image
  suffix = strrep(quantity, 'spheres', '');
  dipole = false;
  if ~isempty(strfind(suffix, 'dipole'))
    dipole = true;
    suffix = strrep(suffix, 'dipole', '');
  end
  energies = str2num(strtrim(suffix));
  if isempty(energies) % default if no list of energies [eV] given after the 'spheres ' prefix
    energies = [2 9 16 23 30];
  end
  if length(energies) <= 2
    n = 20000; % high statistics for few rings (will be doubled if dipole)
  else
    n = 10000; % (still quite a lot) (will be doubled if dipole)
  end
  % Isotropic or dipole-transition like (density proportional to y^2), random directions, no grouping for VMI resolution calculation:
  variables = {'_number_of_directions',iif(dipole,2,1), 'particles_per_direction',n, '_isotropic',1, 'trajectories',true};
  mass = 0; % electron, will get negative charge
  
  variables{end+1} = 'debug'; variables{end+1} = true; % option to keep temporary files for loading to re-plot later (not expecting to run "spheres" often)

else
  % No mass, just electron VMI:
  
  % Customization of which electron energies
  if ~isempty(strfind(quantity, ' h3en')) 
    % Only the higher three energies among the 5 defined.
    lo = [];
  else
    lo = 1;
  end
  
  if ~isempty(strfind(quantity, 'y30_')) % high-energy electron VMI simulation
    energies = [2*lo 9*lo 16 23 30]; %[eV] approximately the same ratios as the low-energy case, but scaled up to 30 eV max (40 eV photon energy - 10 eV ionization energy), 1.82 rounded to 2 which shouold simplify somewhat
  elseif ~isempty(strfind(quantity, 'y42_')) % append also high-energy electron VMI simulation
    energies = [2*lo, 9*lo, 16, 23, 30] *42/30; % [eV] scale to higher energy
    
  elseif ~isempty(strfind(quantity, 'y20_')) % scale the y30 to intermediate energy range, 10-20 eV probably more relevant than 30eV.
    %energies = [1*lo,  4.5*lo,  8,   11.5, 15  ]; % [eV] 3.5*2 eV spacing. Skip intermediate energies to speed up this combined ion + electron optimization.
    energies = [1.5*lo, 6*lo  , 10.5, 15,   19.5]; % [eV] 4.5*2 eV spacing. Skip intermediate energies to speed up this combined ion + electron optimization.
  elseif ~isempty(strfind(quantity, 'y10_')) % scale the y30 to low-energy range, 0-10 eV
    energies = [0.67*lo, 3*lo  , 5.33, 7.67,   10.0]; % [eV] 
  else % original low energies
    if isempty(strfind(quantity, 'y_')) 
      warning('Treating unknown quantity %s as electrons from 0.4 to 6.6 eV.', quantity);
    end
    energies = [0.4*lo, 1.95*lo, 3.5, 5.05, 6.6]; % [eV] 1.55 eV spacing
  end
  
  mass = 0; % electron, will get negative charge
%  variables = {'_number_of_directions',8, 'particles_per_direction',20}; % seems to work well, could possibly increase number of particles for fine-tuning
%  variables = {'_number_of_directions',8, 'particles_per_direction',50}; % incresed 2015-05-28 11:51 for increased reliability
  if quasiisotropic_for_VMI
    variables = {'_number_of_directions',8, 'particles_per_direction',50, '_isotropic',3}; % From 2015-12-30 can use new 3D-isotropic VMI-style, 
    % since 50 instead of 38 particles per direction, the above is not exactly equivalent to what was run in the mass+ combination.
  else
    variables = {'_number_of_directions',8, 'particles_per_direction',50, '_isotropic',2}; % but when eVMI without mass perhaps some interest to cross-validate using old _isotropic=2
  end
end


% variables{end+1} = 'debug'; variables{end+1} = true; % option to keep temporary files for loading to re-plot later (not expecting to run "spheres" often)


% Allow suffix on the form ' vy340' to set initial_vy offset speed [m/s]
match = regexp(quantity, ' vy([0-9]+)', 'tokens');
if ~isempty(match)
  offset_velocity_y = str2double(match{1}{1})
  offset_velocity_y = offset_velocity_y / 1000; % convert from [m/s] to [km/s] used in SIMION & Lua, and change the sign since molecular source will be on top
  variables{end+1} = 'initial_vy'; variables{end+1} = offset_velocity_y;
end
% Allow suffix on the form ' vz340' to set initial_v_vertical offset speed [m/s]
match = regexp(quantity, ' vz([0-9]+)', 'tokens');
if ~isempty(match)
  offset_velocity_z = -str2double(match{1}{1})
  offset_velocity_z = offset_velocity_z / 1000; % convert from [m/s] to [km/s] used in SIMION & Lua, and change the sign since molecular source will be on top
  variables{end+1} = 'initial_v_vertical'; variables{end+1} = offset_velocity_z;
end

% variables{end+1} = 'debug'; variables{end+1} = true; % DEBUG: option to keep temporary files for debugging

%old intent: objective_mode = 0; % return vector, let optmizer take 2-norm (among results for energies) -- not tested, might not be allowed
objective_mode = 1; % average (now 3-norm among results for eVMI-energies, arithmetic between mass and eVMI). Default and probably always used until 2016-01-11.
% Now objective_mode > 1 instead means that pareto frontier search is used, output should have one column per objective. 
% See the function comment in optimize_potentials_objective.m for variants.

save_directory = '.\output\';  % Use a path relative to current directory
% save_directory = 'C:\Users\<USERNAME>\Documents\simulation output\'; % TODO: To be customized (optionally) to relevant path for each user (

%% Handle input arguments

if nargin < 5 || isempty(potential_update_matrix)
  potential_update_matrix = eye(length(potentials_min));
end
if nargin < 6
  absolute_index = [];
end
if nargin < 7
  quantity = 'y_rel mix'; % Optimize for arithmetic mean of y, y2 and y2_sin_weighted relative errors
%   quantity = 'y_rel mix pow'; % Optimize for arithmetic mean of y, y2 and y2_sin_weighted resolving powers
  % quantity = 'y2_rel';
  % quantity = 'y2_rel sin';
end
if nargin < 8
  if isempty(potentials_initial)
    algorithm = 'genetic';
  else
    algorithm = 'active-set'; % default, faster than interior-point but often not quite as good optimum
  end
end
if nargin < 9
  batch_name = ''; % to not write to log file
end
if isempty(algorithm) && (isempty(batch_name) || strcmp(batch_name, 'just evaluate initial simulation'))
  batch_name = 'just evaluate initial simulation';
  algorithm = 'active-set'; % need to put some algorithm to not raise error, value doesn't matter when just evaluating initial
end


if isempty(potentials_initial)
  % Arithmetic mean
  potentials_initial = (potentials_min + potentials_max) / 2;
  
%   % Use geometric mean for non-fixed potentials that have nonzero min & max of the same sign
%   which = potentials_min ~= potentials_max & potentials_min .* potentials_max > 0;
%   % Set initial guess to geometric mean of min and max bounds:
%   potentials_initial(which) = sign(potentials_min(which)) .* sqrt(potentials_min(which) .* potentials_max(which));
end

if isempty(absolute_index)
  relative = false; % directly use the potentials as parameters
else
  relative = true; % probably a more efficient parametrization: One parameter is the A-repeller potential (coefficient wrt. initial guess) and other parameters are corefficients wrt. A-repeller. Thus rescaling lens for higher/lower kinetic energy requires only one parameter to modified.
end



%% Start the work. Not much user-parameters below here.
global optpot_counter
optpot_counter = 0; % reset the global counter used by optimize_potentials_objective for unique filenames
clear optpot_counter;
%assignin('base', 'optpot_counterB', 0); %DEBUG

global optpot_thread_id; optpot_thread_id = []; % clear the calling thread id, because some optimizers makes first evaulation in this thread
global cached_recent; cached_recent = [];

if strcmp(batch_name, 'just evaluate initial simulation')
  % Won't use cached data when returning flys structure to show mass resolution, might as well take it as a moment to ensure optimization cache is cleared
  
  % (Takes several seconds to stop and restart, so skipping it now): matlabpool close force local 
  global cached_recent; cached_recent = []; % quick way to clear the cache used in optimize_potentials_objective.m

else
%elseif false % DEBUG disabled parallellism

  % Allow parallelism (not sure optimizer will use it)

%   try  
%   % First close any old threads, to clear their caches. Kind of slow thouguh, so another way to tell each thread to clear the cache would be useful.
%   % Caller needs to explicitly close the matlatpool first to clear cache. Is done at the top of optimize_potentials_loop.m
%   % If this closing is not done, then cache will persist between optimizations which is bad if external things like SIMION workbench (with same name) is modified.
% %     matlabpool('close');
%      matlabpool close force local
%   catch e
%   end
  
  old_pool = gcp('nocreate'); % in earlier Matlab version: old_size = matlabpool('size');
  if threads == 1 && isempty(old_pool)
    % Don't start any parallelism, not desired
  elseif old_pool.NumWorkers ~= threads
    % Need to adjust size of thread pool
    try
      delete(gcp()); % matlabpool('close');
    catch e
    end
    try
      if threads >= 2
        parpool(threads);
      end
    catch e
      e
    end
  end
end


i = strfind(algorithm, ' ');
if ~isempty(i)
  algorithm_options = algorithm(i(1):end); % include preceding space, so each option word has a preceding space for string matching later
  algorithm = algorithm(1:i(1)-1);
else
  algorithm_options = '';
end
genetic = false; % Genetic algorithm for optimization? If not: gradient-based optimization from inital guess

switch algorithm
  case 'interior-point'
    options = optimoptions('fmincon', 'algorithm', 'interior-point'); % since documentation says this will be default in future, I assume it is generally good
  case 'interior-point cg'
    options = optimoptions('fmincon', 'algorithm', 'interior-point'); % since documentation says this will be default in future, I assume it is generally good
    options.SubproblemAlgorithm = 'cg'; % not sure what it means, but the number of function evaluation drops slightly (1-2%) for interior-point
  case 'active-set'
    options = optimoptions('fmincon', 'algorithm', 'active-set'); % faster than other (about half as many function evaluations) but not always quite as good optimum. Still OK.
  case 'sqp'
    options = optimoptions('fmincon', 'algorithm', 'sqp'); % as slow as interior-point but fewer options, don't use. Well, some recommend it when starting from a good initial point.
  case 'trust-region-reflective'
    options = optimoptions('fmincon', 'algorithm', 'trust-region-reflective'); % current Matlab default, but not availble when some potential held fixed. Will hardly be useful.

  case 'genetic'
    options = gaoptimset('ga');
    genetic = true; % genetic optimization (not requiring initial guess), stochastic result

%   options = psoptimset('patternsearch'); % nearly as fast as active-set for non-relative (not good for relative), but not using gradient (making simple steps in grid, and picking best for next iteration's start)
%   options = saoptimset('simulannealbnd'); % bad
  case {'pareto', 'paretoI', 'paretoIE'}
    % Search for Pareto frontier rather than an individual configuration of the potentials.
    options = gaoptimset('gamultiobj');
    genetic = true; % pareto is a variant of genetic optimization
    % The Pareto frontier contains different ways of (optimally) compromizing between multiple objective, while maintaining spread to not just favour one objective.
    switch algorithm
      case 'paretoI'
        objective_mode = 3; % separate objectives for ion masses (if quantity contains any mass) & at most one for electronVMI
      %Currently not supported: case 'paretoE' separate for electronVMI only)
      case 'paretoIE'
        objective_mode = 4; % separate objectives for each ion mass (if any) and each electron energy (if any)
      otherwise % Default 'pareto'
        objective_mode = 2; % At most two objectives: average of ion masses (if any) & average of VMI for electron energies (if any)
    end
    
  otherwise
    error('Unknown algorithm %s', algorithm);
end


% Definition of SIMION electrode indices for the workbench used (as in Lua workbench script)
if all(size(potentials_initial) >= 2)
  % Matrix given
  if ~genetic
    error('Multiple rows of initial potentials only allowed for genetic algoritihms, not for gradient-based.');
  end
  all_potentials_initial = potentials_initial;
  potentials_initial = mean(potentials_initial); % row-vector of averages used where a single array of potentials needed.
else % row or column vector
  potentials_initial = potentials_initial(:)'; % ensure row-vector
  all_potentials_initial = potentials_initial;
end

potentials_min = potentials_min(:)';
% Take NaN limits to mean "fix at initial guess". Using mean of initial guesses in case multiple given (thus error/outside limits if gusses differ in a "fixed" potential)
potentials_min(isnan(potentials_min)) = mean(potentials_initial(:,isnan(potentials_min)),1);
potentials_max(isnan(potentials_max)) = mean(potentials_initial(:,isnan(potentials_max)),1);

dependent = any(potential_update_matrix~=eye(size(potentials_initial,2)),2)'; % which potentials will be set to a linear combination of other?
tmp = dependent & potentials_min ~= potentials_max; % any conflict?
if any(tmp)
  error('The following potentials will be altered by the potential_update_matrix but are also marked for optimization (within a range).\n Indices in potential array: %s', num2str(find(tmp)))
end


if ~genetic
  
  % Stopping criteria
  if ~isempty(strfind(algorithm_options, ' quick-debug'))
    options.MaxIter = 1; options.MaxFunEvals = 12; %for debugging only
  elseif ~isempty(strfind(algorithm_options, ' quick'))
    options.MaxIter = 15; options.MaxFunEvals = 50; %for rough start
  elseif ~isempty(strfind(algorithm_options, ' fine'))
    options.MaxIter = 25; options.MaxFunEvals = 200; % probably more than necessary, better to try other intial guesses than optimize this far
  else % default
    options.MaxIter = 20; options.MaxFunEvals = 130; %for OK/good quality
  end
  if ~relative
%     options.TolX = 0.00001 / max(max(abs(potentials_initial))); % (relative!) must be quite small, in order to not stop before first step
    options.TolX = 0.000005 / max(max(abs(potentials_initial))); % (relative!) decreasing to try harder also when very good initial guess (2015-06-17)
  else
%    options.TolX = 0.001; % (relative!) works better when using relative parametrization
   options.TolX = 0.0005; % (relative!) decreasing to try harder also when very good initial guess (1V/2000V)
  end
  options.PlotFcns = {@optimplotx, @optimplotfval,@optimplotstepsize};
    
else % genetic
  
  % For genetic
  if ~isempty(strfind(algorithm_options, ' quick-debug'))
    options.Generations = 2; options.PopulationSize = 6; options.EliteCount = 2; % too small to be good
    options.StallGenLimit = 3;
  elseif ~isempty(strfind(algorithm_options, ' quick'))
    %options.Generations = 10; options.PopulationSize = 8; options.EliteCount = 3; % probably close to minimum size, can get stuck
    options.Generations = 7; options.PopulationSize = 8; options.EliteCount = 3; % probably close to minimum size, can get stuck
    options.StallGenLimit = 5; % to stop earlier than default 50 (or max 10) if progress is lacking for 5 generations
  elseif ~isempty(strfind(algorithm_options, ' fine'))
    options.Generations = 15; options.PopulationSize = 16; options.EliteCount = 4; % when a bit more time can be used
    options.StallGenLimit = 15; % to not stop earlier just because (random) progress is lacking
  elseif ~isempty(strfind(algorithm_options, ' many')) % many individuals, moderate number of generations. Guessing useful for Pareto search (spending more time but not needing many restarts to get diversity of results)
    options.Generations = 6; options.PopulationSize = 30;
    options.EliteCount = 9; options.ParetoFraction = 9/30; % Just below default ParetoFraction=0.35. Probably OK as it is only an upper limit, not certain that this many "on front" will be found.
    options.EliteCount = 7; options.ParetoFraction = 7/30; % EliteCount not used by Pareto, approach that effect slightly by lowering ParetoFraction from default 0.35
    options.EliteCount = 7; options.ParetoFraction = 9/30; % 2016-04-06 Raising ParetoFraction again, since finding that often only one of the 7 has a non-miss for third objective (42eV electron), and adding fourth objective for pot#10.
    options.StallGenLimit = 6; % may stop a bit earlier if no progress
  elseif ~isempty(strfind(algorithm_options, ' qmore')) % more individuals and fewer generations than 'qmany', hopefully better when giving many good initial guesses and only small improvement expected
    %options.Generations = 6; options.PopulationSize = 17; options.EliteCount = 5; % product 6*17 is smaller than qmany's 9*12 so possibly a bit faster
    options.Generations = 6; options.PopulationSize = 17; options.EliteCount = 4; % 2016-01-06 trying with fewer elite, for diversity to survive longer
  elseif ~isempty(strfind(algorithm_options, ' qmany7')) % intermediate between quick and 'fine', not many generations but more individuals
    %options.Generations = 7; options.PopulationSize = 12; options.EliteCount = 4; % speed up more compared to 'qmany', 7 generations usually enough
    options.Generations = 7; options.PopulationSize = 12; options.EliteCount = 3; % 2016-01-06 trying with fewer elite, for diversity to survive longer
    options.StallGenLimit = 10; % to not stop earlier just because (random) progress is lacking
  elseif ~isempty(strfind(algorithm_options, ' qmany')) % intermediate between quick and 'fine', more emphasis on individuals than generations. Preferred choice.
    options.Generations = 9; options.PopulationSize = 12; options.EliteCount = 4;
    options.StallGenLimit = 10; % to not stop earlier just because (random) progress is lacking
  else % default
    options.Generations = 10; options.PopulationSize = 10; options.EliteCount = 3; % OK compromize
    options.StallGenLimit = 7; % to stop earlier than default 50 (or max 10) if progress is lacking for 7 generations
  end
  options.PopInitRange = [0.4; 1.1]; % (for relative parametrization) (will be expanded with identical columns per parameter)
  %options.FitnessLimit = 1E-3; % stops if this low fitness would be reached
  options.FitnessLimit = 5E-4; % stops if this low fitness would be reached
  
  % NOTE: with FitnessScalingFcn = @fitscalingrank (default), the direct fitness values are not used. Instead results are sorted and 1/sqrt(rank) is used.
  % This means my high values for particle-misses don't make individuals as unlikely as thought, but is probably OK (equalizing the score for all the bad individuals).
  
%   options.SelectionFcn = @selectionstochunif; % probability ~ 1/scaledfitness, unlikely to miss individual with very good 1/scaledfitness. Default. Used always until 2016-01-06.
%   options.SelectionFcn = @selectionroulette; % probability ~ 1/scaledfitness, seems like most fair random with large population but could miss several of the best by chance
%   options.SelectionFcn = @selectiontournament; % based on small-group comparisons rather than probability ~ 1/scaledfitness, e.g. size=4 means worst 3 individuals will never be chosen but any other may be (with higher chance for the good ones)
  
%   options.CreationFcn = @gacreationuniform; % Matlab default
%   options.CreationFcn = @gacreationlinearfeasible; %DEDBUG
%   options.CreationFcn = 'gacreationnonlinearfeasible'; %DEDBUG
  
  if ~isempty(strfind(algorithm_options, ' mut:uniform'))
    options.MutationFcn = {@mutationuniform, 0.04}; % uniform instead of gaussian, to not violate bounds so easily. Raise the mutation rate from default 0.01 in uniform.
  elseif ~isempty(strfind(algorithm_options, ' mut:adapt'))
  	options.MutationFcn = @mutationadaptfeasible; % didn't seem helpful (but only tested with two variable potentials)
  else % default
    %options.MutationFcn = {@mutationgaussian, 0.4, 0.95}; % Good although children can be generated outside bounds. The current parameters give lower initial std (since values near bounds are unlikely to be good). Also a lower the 'shrink' from default 1, to retain some mutation spread also in last iteration
    options.MutationFcn = {@mutationgaussian, 0.5, 0.85}; % Reduced shrink further to get less stuck with multiple identical individuals, and incresed initial spread a bit.
  end
  options.CrossoverFraction = 0.67; % reduce crossover a bit (default 0.8), in favour of more mutations to keep diversity high
  
  if ~isempty(strfind(algorithm_options, ' cro:scatt'))
    options.CrossoverFcn = @crossoverscattered; % OK, Matlab default
  elseif ~isempty(strfind(algorithm_options, ' cro:heur'))
    options.CrossoverFcn = @crossoverheuristic; % seems this gets "stuck" easier, seemingly not as good. (In theory the step to and beyond better parent seems like it should help avoid getting stuck, though.)
  else % default
    options.CrossoverFcn = {@crossoverintermediate, ones(1,size(potentials_initial,2))}; % OK, preferred. Can still get stuck if population spread becomes too small.
  end
  
  if objective_mode < 2
    options.PlotFcns = {@gaplotrange, @gaplotdistance, @gaplotgenealogy,};
  else % Pareto
    if ~isempty(strfind(algorithm_options, ' genotype'))
      % By default the distance measurement (controlling diversity) is based on phenotype, i.e. on the spread of objectives. (Perhaps that means it thinks keeping some of the really bad ones with high badness due to misses is good?!)
      % Change it to genotype, i.e. based on the spread of parameters for the potentials.
      options = gaoptimset(options, 'DistanceMeasureFcn', {@distancecrowding,'genotype'});
    end
    
    % DEBUG which plots are useful?
    options.PlotFcns = {@gaplotrankhist, @gaplotdistance, @gaplotpareto, @gaplotparetodistance, @gaplotspread, @gaplotgenealogy};
    % IMPROVEMENT  if @gaplotpareto only plots first two, so could do variant where both masses combined to one score vs eVMI score.
    %  Could make it as a custom plot function http://uk.mathworks.com/help/gads/genetic-algorithm-options.html#f15210
  end

end

% For patternsearch (seems OK for very rough trend, but slow to reach good optimum and perhaps too local)
% options.Cache = 'on'; options.CacheTol = 1E-3;
% options.ScaleMesh = 'off';
% %options.TolX = 0.5; % options.InitialMeshSize = 50; %[V] % absolute
% options.InitialMeshSize = 0.15; % relative
% non-relative % Optimized potentials: -1449.2 	-1449.2 	 -958.4 	    0.0 	    0.0 	    0.0	: 0.012163 = 1/82.219, func count 55
% patttern initial 0.5, scale, cache 1E-4    Optimized potentials: -534.37  -534.37  -395.44     0.00     0.00     0.00	: 0.151114 = 1/6.618 (in 54 evaluations)
% patttern initial 0.1, no scale, cache 1E-4 Optimized potentials: -532.50  -532.50  -394.05     0.00     0.00     0.00	: 0.150943 = 1/6.625 (in 54 evaluations)

if ~isempty(strfind(quantity, '_rel')) || ~isempty(strfind(quantity, 'mass'))
  % a relative error (inverse resolving power), good values are < 0.03 (3%, resolving power 33)
  % note: if gradient evaulation over options.DiffMinChange step doesn't cause larger change than TolFun, optimizer stops
  if ~isempty(strfind(algorithm_options, ' quick'))
  	options.TolFun = 2E-4; % OK, somewhat rough
  else
%    options.TolFun = 5E-5; % when aiming higher  
    options.TolFun = 5E-6; % when aiming higher. Lowered 2015-06-17 since for (good) mass resolution, the non-genetic methods never make more than one step (mostly due to local minimum, perhaps more about step length)
  end
else
  % optimizing an absolute standard deviation ([mm], [mm^2] or [ns]) --- IMPROVEMENT: not implemented/tested yet
  if ~isempty(strfind(algorithm_options, ' quick'))
    options.TolFun = 0.05;
  else
    options.TolFun = 0.01;
  end
end

% Other options
options.UseParallel = 'always';
options.Display = 'iter';
  
if any(strcmp(variables, 'potentials'))
  error('No ''potentials'' may be set in shared array of variables when optimizing.');
end
PUMT = transpose(potential_update_matrix); % since potentials is a row-vector rather than column-vector

if relative
  % Try to reformulate parameters to have a scaling factor, for easier adaptation of all potentials (scaling preserves the focal length)
  %absolute_index = B_FIRST; % let other potentials be expressed as coefficients of the A-repeller potential
  potentials_abs = mean(all_potentials_initial,1); % to have a scalar reference potential (even when several initial potentials for genetic population given), take mean among rows
  if any(potentials_abs(:,absolute_index) == 0)
    error('May not use zero as initial guess for the scaling potential at index %d.', absolute_index);
  end
  options.TypicalX = max(mean(abs(all_potentials_initial),1), 100) ./ abs(potentials_abs(absolute_index)); % provide info about typical values (magnitude)
  potentials_initial = potentials_initial ./ potentials_abs(absolute_index); % the initial guess will always be 1.0 for the absolute_index potential coefficient
  
  % Logic error discovered 2015-12-18 for relative parameters when multiple initial guesses:
  % Old only correct for single initial guess (not handling the fact that potentials_abs is average of all initial guesses):
  %   all_potentials_initial = all_potentials_initial ./ potentials_abs(absolute_index); % the initial guess will always be 1.0 for the absolute_index potential coefficient
  % Correction:
  tmp = all_potentials_initial ./ repmat(all_potentials_initial(:,absolute_index), 1,size(potentials_initial,2)); % get the relative coordinates within each initial guess, here the columns for absolute_index is always 1.0.
  tmp(:,absolute_index) = all_potentials_initial(:,absolute_index) ./ potentials_abs(absolute_index); % now set each initial guess relative to the single absolute reference
  all_potentials_initial = tmp;
  % The old code resulted in each initial guess being scaled from the average by an extra factor equal to its initial_potential(absolute_index)/potentials_abs(absolute_index)
  % which is not a big problem when all the intials have similar values for the abslute index.
  % This bug explains why I this fall stopped using the relative mode, while when today trying it (values for all initials within 1%) it worked well again.


  potentials_abs_min = potentials_min; potentials_abs_max = potentials_max;
  %if narrow
  % Translate absolute min & max in to relative, in a narrow way (only using ratio in case the absolute_index coefficient is 1.0
  %   potentials_min = potentials_min ./ potentials_abs(absolute_index);
  %   potentials_max = potentials_max ./ potentials_abs(absolute_index);
  % else % Alt. "worst-case" ratio using min/max of the absolute_index coefficient. This means other parameters may be allowed outside the original non-relative ranges
%     potentials_min = potentials_min ./ potentials_abs_max(absolute_index);
%     potentials_max = potentials_max ./ potentials_abs_min(absolute_index);
  % else % Alt. heuristic weighting to not give as extreme range as "worst case" (to benefit from the relative scaling aspect also in the limits, e.g. important for initial random population of genetic algorithm)
    potentials_min = potentials_min ./ [potentials_abs_max(absolute_index)*2 + potentials_abs(absolute_index)*1]*3;
    potentials_max = potentials_max ./ [potentials_abs_min(absolute_index)*2 + potentials_abs(absolute_index)*1]*3;
  % end
  
  if potentials_abs(absolute_index) < 0
    % min & max are swapped if divided by negative number
    if potentials_abs_max(absolute_index) > 0 || potentials_abs_min(absolute_index) > 0
      % DEBUG whether min&max with differnet sign can be translated correctly
      error('Translation of min and max bounds to relative parameters is not supported when the main potential''s (index %d) bounds have different signs, i.e. negative as well as positive allowed.', absolute_index);
    end
    tmp = potentials_max;
    potentials_max = potentials_min;
    potentials_min = tmp;
  end
  
  % A simple parametrization that is just rescaled in magnitude, without a common scaling parameter, would be:
  %pot = @(param) potentials_abs(absolute_index) * param;
  % With global scaling parameter, and relative coefficients for other potentials:
  powers = ones(1,size(potentials_initial,2));
  powers(absolute_index) = 0;
  pot = @(param) potentials_abs(absolute_index) * param(absolute_index) * param.^powers;
  
  model = @(param) optimize_potentials_objective(workbench_name, mass,energies, source_point_distribution, ...
                     [variables, {'potentials', pot(param)*PUMT}], quantity, objective_mode);
  if ~genetic
    %options.DiffMinChange = 1E-4; % The smallest step to take when evaluating gradient. If simulation is noisy or rounds potentials, this should not be too small.
    options.DiffMinChange = 2E-4; % The smallest step to take when evaluating gradient. If simulation is noisy or rounds potentials, this should not be too small.

    % Try to enforce larger (some volt) step sizes than default, since rounding errors (noise) make a comparison of SIMION result for nearly identical parameters uncertain:
    %options.DiffMinChange = min(min(1E-2, 2.0 ./ abs(potentials_abs(potentials_min<potentials_max)))); % ~ what corresponds to 2 V (or 1% of typical value when value is near zero)
%     options.DiffMinChange = min(min(5E-3, 1.0 ./ abs(potentials_abs(potentials_min<potentials_max)))); % ~ what corresponds to 1 V (or 0.5% of typical value when value is near zero). Was used until (including) 2015-05-21.
    tmp = min(5E-3, 2.0 ./ abs(potentials_abs(potentials_min<potentials_max))); options.DiffMinChange = mean([mean(tmp) min(tmp)]); % don't let parameter with highest |potentials_abs| fully decide, let average value influence too.
    
%    options.DiffMaxChange = max(min(0.25, 25 ./ abs(potentials_abs(potentials_min<potentials_max)))); % ~ what corresponds to 25 V (or 25% of typical value when value is near zero)
    options.DiffMaxChange = max(min(0.25, 30 ./ abs(potentials_abs(potentials_min<potentials_max)))); % ~ what corresponds to 3 V (or 25% of typical value when value is near zero)
    % options.FinDiffRelStep = 5 / mean(abs(potentials_initial)); % ? scale all derivative steps up by this factor
    options.FinDiffRelStep = min(1E-2, 1.0 ./ abs(potentials_abs)); % ~ what corresponds to 1 V (or 1% when value is near zero)
    % options.InitBarrierParam = 0.5; % ? -- Leaving at default 0.1
  end
  
else % Simply use the potentials as parameters
  potentials_abs = potentials_initial; % update potentials_abs anyway
  model = @(pot) optimize_potentials_objective(workbench_name, mass,energies, source_point_distribution, ...
                     [variables, {'potentials', pot*PUMT}], quantity, objective_mode);
  if ~genetic
    options.DiffMinChange = 2; %0.5; % The smallest step to take when evaluating gradient. If simulation is noisy or rounds potentials, this should not be too small.
    options.DiffMaxChange = 25;
  end
  %(original, could set negative TypicalX) tmp=potentials_initial; tmp(tmp==0) = 50; options.TypicalX = tmp; % provide info about typical values (same as initial guess, unless guess potential is zero)
  options.TypicalX = max(mean(abs(potentials_initial),1), 50); % provide info about typical values (magnitude)
end

if ~genetic && ~isempty(strfind(algorithm_options, ' wide'))
  % Aim for exploring wider landscape or managing noisy simulation, avoid small steps
  options.DiffMinChange  = options.DiffMinChange  * 6;
  options.DiffMaxChange  = options.DiffMaxChange  * 12;
  options.FinDiffRelStep = options.FinDiffRelStep * 12;
  options.TolX = options.TolX / 5;  % don't know if relevant but should reduce risk of stopping at initial point
  options.TolFun = options.TolFun / 10;  % don't know if relevant but should reduce risk of stopping at initial point
end

if strcmp(batch_name, 'just evaluate initial simulation')
  % A simple way to call runsim, benefiting from the potential_update_matrix used in this function.
  if size(potentials_initial,1) > 1
    error('Running/showing multiple simulations not handled by optimize_potentials. See caller (optimize_potentials_loop.m).');
  end
  potentials = potentials_initial * PUMT;
  [lowest_cost, costs_before_averaging, flys] = model(potentials_initial);
  lowest_cost = costs_before_averaging; % may return array in the 'just evaluate initial simulation' case
  output = flys; % give this output argument a different meaning
  disp('NOTE: Only evaluating for initial potentials.');
  return;
end

if isempty(strfind(algorithm_options, ' keep_constrained_dimensions'))
  % (default) Shrink the problem to only keep non-fixed potential parameters.
  % It seems fmincon otherwise evaluates gradient also in the non-variable dimensions, which wastes time.
  small = true;
  
  optimized_indices = find(potentials_min < potentials_max);
  if size(potentials_initial,1) > 1
    error('Running/showing multiple simulations currently not handled when keep_constrained_dimensions.');
  end
  small_model = @(param) model(inplace(optimized_indices, param, potentials_initial));
  if length(options.TypicalX) > 1
    options.TypicalX = options.TypicalX(optimized_indices);
  end
  if isprop(options, 'FinDiffRelStep') && length(options.FinDiffRelStep) > 1
    options.FinDiffRelStep = options.FinDiffRelStep(optimized_indices);
  end
  if isfield(options, 'CrossoverFcn') && iscell(options.CrossoverFcn) && strcmp(char(options.CrossoverFcn{1}),'crossoverintermediate')
    options.CrossoverFcn{2} = options.CrossoverFcn{2}(optimized_indices);
  end
 
else % original version, keeping also the fixed parameters
  small = false;
  small_model = model;
  optimized_indices = 1:size(potentials_initial,2);
end

if threads < 2
  options.UseParallel = false;
end

tic
if ~genetic
  %[final_parameters,lowest_cost,exitflag,output,lambda] = fmincon(model, potentials_initial, [],[], [],[], potentials_min, potentials_max, [], options);
  [final_parameters,lowest_cost,exitflag,output] = fmincon(small_model, potentials_initial(1,optimized_indices), [],[], [],[], potentials_min(optimized_indices), potentials_max(optimized_indices), [], options);
else

  % It is OK (at least for the @gacreationuniform) to provide a PARTIAL initial population, e.g. containing the initial guess!
  % This should avoid the problem that with a small fully random population, most initial guesses get no particles reaching the detector.
  if relative
    tmp = inplace(absolute_index, 1, zeros(1,size(potentials_initial,2))); % only scale the common parameter (to retain presumably acceptable lens shape of initial guess)
  else
    tmp = ones(1,size(potentials_initial,2)); % scale all parameters (achievning the same result as in relative case)
  end
  
  % Generate initial population
  options.InitialPopulation = all_potentials_initial(:,optimized_indices); % Use the given inital potential array, or set of potential arrays
  
  if size(all_potentials_initial,1) ~= 1 && size(options.InitialPopulation,1)+2 < options.PopulationSize % (ensure at least two remains to be generated by randomization below or options.CreationFcn)
    options.InitialPopulation(end+1,:) = potentials_initial(optimized_indices); % add also the average
  end
  
  if size(options.InitialPopulation,1)+iif(relative,1,2) < options.PopulationSize % (ensure at least one (if relative) or two (otherwise) remains to be generated by options.CreationFcn)
    %options.InitialPopulation(end+1,:) = potentials_initial(optimized_indices) .* 1.5.^tmp(optimized_indices); % scale all potentials of initial guess by 1.5
    options.InitialPopulation(end+1,:) = potentials_initial(optimized_indices) .* 1.2.^tmp(optimized_indices); % scale all potentials of initial guess by 1.2
  end
  if relative && size(options.InitialPopulation,1)+2 < options.PopulationSize % (ensure at least two remains to be generated by randomization below or options.CreationFcn)
   %options.InitialPopulation(end+1,:) = potentials_initial(optimized_indices) .* 0.5.^tmp(optimized_indices); % halve all potentials of initial guess
   options.InitialPopulation(end+1,:) = potentials_initial(optimized_indices) .* 0.8.^tmp(optimized_indices); % scale all potentials of initial guess by 0.8
  end
  
  for i = 1:3
    if size(options.InitialPopulation,1)+1 < options.PopulationSize % (ensure at least one remains to be generated by options.CreationFcn)
      % Generate one individual using a triangular distribution, centered in the middle of the allowed range
      % then take weighted average of it and the initial guess.
      options.InitialPopulation(end+1,:) = (potentials_min(optimized_indices) ...
        + (mean(rand(2,length(optimized_indices)),1) .* (potentials_max(optimized_indices)-potentials_min(optimized_indices)))) * 0.5 ...
        + potentials_initial(optimized_indices) * 0.5;
    end
  end
  
  if objective_mode < 2 % not Pareto search, but normal to get a single "optimal" configuration
    %[final_parameters,lowest_cost,exitflag,output] = ga(model, length(potentials_initial), [],[], [],[], potentials_min, potentials_max, [], options);
    if relative
      [final_parameters,lowest_cost,exitflag,output] = ga(small_model, length(optimized_indices), [],[], [],[], potentials_min(optimized_indices), potentials_max(optimized_indices), [], options);
    else
      % 2016-01-06 Learned how to specify that parameters are integer! (Give their indices)
      %disp('DEBUG: trying with integer [V] potentials!'); % DEBUG verify that it works. IMPROVEMENT Try to compare if it is better / less stuck after many generations?
      [final_parameters,lowest_cost,exitflag,output] = ga(small_model, length(optimized_indices), [],[], [],[], potentials_min(optimized_indices), potentials_max(optimized_indices), [], 1:length(optimized_indices), options);
    end
    
  else
    % Search for Pareto frontier (different trade-offs/compromises between multiple objectives) instead of single result
    % For this, the objective function should return multiple columns (one per objective).
    
    disp(sprintf('Searching for Pareto frontier, objective_mode=%d.', objective_mode));
%     global final_parameters lowest_cost exitflag output % DEBUG
    [final_parameters,lowest_cost,exitflag,output] = gamultiobj(small_model, length(optimized_indices), [],[], [],[], potentials_min(optimized_indices), potentials_max(optimized_indices), options);
    % Here the result varaibles differ from the usual (single-individual result):
    % * final_parameters is a matrix with one row per "individual"/point on the Pareto frontier, and one column per potential-parameter.
    % * lowest_cost      is a matrix with one row per -''-, and one column for objective.
    % * output.averagedistance may say something about how well the points are spread along the frontier
    % * output.spread may do the same or more about how much the frontier moved by the last generation
    
    disp('Pareto frontier individuals');
    %num2str([final_parameters lowest_cost]
    disp(num2str(final_parameters, '%.4g,'));
    disp(num2str(lowest_cost, '%.4g '));
    % DEBUG...
  
    % TODO IMPROVEMENT: PlotFcn @gaplotpareto only plots first two.
    % Make a figure with multiple 2D-subfigures with each combination of two of the objectives, and use different markers for each row pareto point.
    % Could also do a plot variant where mass-objectives are averaged and electron-objectives averaged, to get a 2D plot even if each mass or energy was a separate objective.
    
  end
  
end
elapsed = toc;

%   [final_parameters,lowest_cost,exitflag,output] = patternsearch(model, potentials_initial, [],[], [],[], potentials_min, potentials_max, [], options);
potentials   = NaN(size(final_parameters,1), size(all_potentials_initial,2));
if small
  tmp        = NaN(size(final_parameters,1), size(all_potentials_initial,2)); % can't overwrite in final_parameters (when more than one result of Pareto) since number of column changes
  for i = 1:size(final_parameters,1)
    tmp(i,:) = inplace(optimized_indices, final_parameters(i,:), potentials_initial);
  end
  final_parameters = tmp; clear tmp
end
for i = 1:size(final_parameters,1)
  if relative
    potentials(i,:) = pot(final_parameters(i,:)) * PUMT;
  else
    potentials(i,:) = final_parameters(i,:) * PUMT;
  end
end
if isfield(output, 'funccount')
  output.funcCount = output.funccount;
  elapsed_per_evaluation = elapsed / output.funccount;
else
  elapsed_per_evaluation = NaN;
end

%str = sprintf(' Optimum %s: %s %.6f = 1/%.3f (in %d evaluations)', sprintf(' %8.2f ', potentials), quantity, lowest_cost, 1./lowest_cost, output.funcCount);
%str = sprintf('%s,  %.6f %%=%s= 1/%.3f (in %d evaluations)', sprintf(' %8.2f ', potentials), lowest_cost, quantity, 1./lowest_cost, output.funcCount);
% 1V precision is enough for the potentials:
if objective_mode < 2 % not Pareto search, but normal to get a single "optimal" configuration
  str = sprintf('%s,  %.6f %%=%s= 1/%.3f (in %d evaluations, %.1f min, %.0f s/eval.). Printing "resolving power" = 1/cost here:', sprintf('%5.0f ', potentials), lowest_cost, quantity, 1./lowest_cost, output.funcCount, elapsed/60, elapsed_per_evaluation);
else
  str = sprintf('%% Found %d points to approximate Pareto front (in %d evaluations, %.1f min, %.0f s/eval.). Printing "resolving power" = 1/cost here:', size(potentials,1), output.funcCount, elapsed/60, elapsed_per_evaluation);
  for i = 1:size(potentials,1)
    str = sprintf('%s\n%s%%%s =1/%s. ', str, sprintf('%5.0f ', potentials(i,:)), num2str(1./lowest_cost(i,:), '%#.4g '), quantity);
  end
end

disp(str)
if ~isempty(batch_name)
  % Write summary to log file

  if ~exist(save_directory, 'file')
    [p,d] = fileparts(save_directory);
    mkdir(p,d);
  end
  f = fopen(fullfile(save_directory, [batch_name '.txt']), 'at');
  
  fprintf(f, '%% Workbench: %s; %s; Optimizing %d/%d potentials %s;', workbench_name, batch_name, ...
      sum(potentials_min < potentials_max), length(potentials), iif(relative, sprintf('relative #%d',absolute_index), 'separately') );
  for i = find(dependent)
    fprintf(f, ' #%d@%s', i, mat2str(potential_update_matrix(i,:),4));
  end
  fprintf(f, '\n');
  
  if genetic
    if iscell(options.MutationFcn)
      mutation = char(options.MutationFcn{1});
    else
      mutation = char(options.MutationFcn);
    end
    if iscell(options.CrossoverFcn)
      crossover = char(options.CrossoverFcn{1});
    else
      crossover = char(options.CrossoverFcn);
    end
    if objective_mode < 2 % not Pareto search, but normal to get a single "optimal" configuration
      fprintf(f, '%% Genetic %d/%d generations, population %d, elite %d, mutation %s, crossover %s\n', output.generations, options.Generations, options.PopulationSize, options.EliteCount, mutation, crossover);
    else
      if ~isempty(strfind(algorithm_options, ' genotype'))
        spread_type = 'genotype';
      else
        spread_type = 'phenotype';
      end
      fprintf(f, '%% Pareto genetic %d/%d generations, population %d, frontier %.1f, spread %s, mutation %s, crossover %s\n', output.generations, options.Generations, options.PopulationSize, options.ParetoFraction*options.PopulationSize, spread_type, mutation, crossover);
    end
    clear mutation crossover
  else
    fprintf(f, '%% Starting point %s, %s, %d/%d iter, TolX %.4g, DiffMin %.3g, SPA %s\n', mat2str(potentials_abs*PUMT,7), options.Algorithm, output.iterations, options.MaxIter, options.TolX, options.DiffMinChange, options.SubproblemAlgorithm);
    %mat2str(options.FinDiffRelStep,3),  output.message
  end
  fprintf(f, '%s\n\n', str);
  fclose(f);
  clear dependent f
end

