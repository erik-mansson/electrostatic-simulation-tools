% Objective function for optimizer in optimize_potentials.
% The main purpose is to evaluate the result into the scalar or vectorial quantity
% that the optimizer will minimize.
% PARAMETERS
%   ... see runsim3
%   quantity          What quantity should be minimized?
%                     'y_std', 'y2_std_mean', 't_std':                               concern the absolute standard deviation ([mm] or [mm^2] or [ns])
%                     'y_rel', 'y2_rel', 'y2_rel sin', 'y_rel mix', 'y_rel mix pow': concern the transverse 1/(resolving power) (dimensionless, y/sy, y^2/sy^2 etc.)
%                     't_rel', 'tot_rel', 'mass':                                    concern other 1/(resolving powers) (dimensionless, t/st, combined, m/range(t))
%                     TODO: implement more of these suggestions.
%                     NOTE: the inverse of resolving powers (i.e. relative error) may be better for the optimizer,
%                     since it more naturally makes the gradient small (small steps) near the optimum.
%
%   objective_mode    Should a single cost or separate costs for multiple objectives be returned?
%                   1:  Arithmetic average between the objectives. When combining mass resolution and electron-VMI,
%                       arithmetic_mean( arithmetic_mean(mass_costs), mean_using_3norm(electron_VMI_costs) ) is used.
%                   2:  Separate objectives for ion mass resolution and electron VMI, i.e. two objectives if both mentioned in the quantity string.
%                       Averaging as usual between the ion masses, and between the electron energies.
%                >= 3:  One objective per mass_cost. And one additional for electron-VMI if included.
%                >= 4:  Separate objectives for everything (each mass_cost, each electron ennergy).
%
%              (old argument "group_avg_norm" had this similar intent but only case 1 was used:
%                     If multiple groups are flewn (e.g. different energies), how to combine their results? 
%                     0: Don't combine, output a vector (which optimizer will optimize 2-norm of).
%                        This gives the optimizer more info than when a scalar is returned, but is not supported by fmincon (optimizer used when some potential is held fixed).
%                     1: Arithmetic average
%                     2: 2-norm (like 0 but performed here, so it is supported when some potential is held fixed)
%               )
% RETURNS
%  badness                 The total cost value (lower is better).
%  costs_before_averaging  Optionally, one or more values, whose average (using objective_mode) is the badness.
%  flys                    Optionally, if a second output argument is taken it retuns simulation structure "flys" from readsim3.
% 
% SEE ALSO
% optimize_potentials.m, optimize_potentials_loop.m, runsim3.m, readsim3.m
% 
function [badness, costs_before_averaging, flys] = optimize_potentials_objective(workbench_name, mass,energies, source_point_distribution, variables, quantity, objective_mode)

if nargin < 6
  quantity = 'y_rel mix';
end
if nargin < 7
  objective_mode = 0;
end
% if nargin < 8
%   store_result = false;
% end

if strcmp(variables{end-1},'potentials')
  % Round to 0.1 V precision, before simulating and checking cache
  variables{end} = round(variables{end}*10)/10;
end

global optpot_counter %opt_pot_cache
if isempty(optpot_counter)
  optpot_counter = 1;
  %opt_pot_cache = {};
else
  optpot_counter = optpot_counter + 1;
end
counter_local = optpot_counter; %clear optpot_counter; % attempt at parallellism

% if evalin('base', '~exist(''optpot_thread_id'',''var'')')
%  assignin('base', 'optpot_thread_id', now()-736125 + round(rand * 1E4)*10 + round(rand * 1E3)*1E-3); % generate an presumably different number for each thread, stored in its base workspace
%   optpot_thread_id = evalin('base','optpot_thread_id');
% end
global optpot_thread_id % assigignin(base) is local to thread just like global, so simpler to use "global" syntax
if isempty(optpot_thread_id)
  optpot_thread_id = now()-736125 + round(rand * 1E4)*10 + round(rand * 1E3)*1E-3; % generate an presumably different number for each thread, stored in its base workspace
  try
    t = getCurrentTask(); optpot_thread_id = optpot_thread_id + t.ID * 10; % Use Matlab builtin (fails in this version??) for getting thread id instead of just the heuristics
  catch e
  end
end
counter_local_raw = counter_local;
counter_local = counter_local + 5*(mod(floor(optpot_thread_id)*10 + floor(mod(optpot_thread_id,1)*1E4) + floor(mod(optpot_thread_id,1E-4)*1E9), 1E4)); % convert id to 5*integer(0 to 9999) and add to the counter

% variables{end+1} = 'debug'; variables{end+1} = true; global fly % DEBUG

% Caching of recent evluations, to avoid repeated SIMION calls e.g. once per generation
% for surviving individuals in in genetic optimization.
%
% When matlabpool with multiple workers used, globals as well as assigning(base) use separate workspaces and the threads do not communicate.
% It means each thread will have its own cache, surviving until next optimization unless matlabpool('close') is called in between,
% but still more efficient (when genetic optimization) than if no caching at all.

% Check if same parameters recently evaluated, using a cace of limited size
global cached_recent
if isempty(cached_recent)
  cached_recent = struct('checksums',{}, 'all_vars',{}, 'badnesses',{});
end
checksum = 0;
all_vars = [{workbench_name, mass, energies, source_point_distribution, quantity, objective_mode}, variables];
for i = 1:length(all_vars)
  if ischar(all_vars{i}) % string
    tmp = double(all_vars{i}) .* (1+mod(double(1):length(all_vars{i}), 37)) / 40430;
  else % number or array of numbers
    tmp = all_vars{i};
    tmp(isnan(tmp)) = mod(checksum, 135+i);
    tmp(2:2:end,:,:) = -0.69867 * tmp(2:2:end,:,:);
  end
  tmp = tmp(:);
  tmp(1:4:end) =  (0.9731+mod(i,3.11)) * tmp(1:4:end);
  tmp(2:4:end) =  (0.8287+mod(i,2.96)) * tmp(2:4:end);
  tmp(3:4:end) = (-1.0291+mod(i,0.60)) * tmp(3:4:end);
  tmp(4:4:end) = (-0.7382+mod(i,4.11)) * tmp(4:4:end);
  tmp = tmp .* ((40 + (1:length(tmp))')/50);
  checksum = checksum + mod(sum(tmp), 64355.61 + 3.531*i) / sum(mod(i,29)*0.1+7+1.3*size(all_vars{i}));
end
%disp(sprintf('Thread %.8f, counter %d (%d). Checksum %.16g', optpot_thread_id, counter_local, counter_local_raw, checksum)) % DEBUG
if nargout == 1 && ~isempty(cached_recent)
  i = find(cached_recent(1).checksums == checksum, 1);
  if ~isempty(i) && isequal(cached_recent(1).all_vars{i}{end}, all_vars{end}) && isequal(cached_recent(1).all_vars{i}, all_vars)
    % all_vars{end} is typically the array of potentials, which varies most
    % Already know result
    badness = cached_recent(1).badnesses(i,:);
    disp(sprintf('Thread %.8f, counter %d (%d). Retrieving %.16g from #%d/%d in thread cache. Fitness:%s', optpot_thread_id, counter_local, counter_local_raw, checksum, i, length(cached_recent(1).checksums), sprintf(' %.6f', cached_recent(1).badnesses(i,:)))); % DEBUG

    % Put this entry last, to keep it counted as recent and not removed later
    cached_recent(1).checksums(end+1,1) = cached_recent(1).checksums(i);
    cached_recent(1).all_vars{end+1,1} = cached_recent(1).all_vars{i};
    cached_recent(1).badnesses(end+1,:) = cached_recent(1).badnesses(i,:);
    % Remove the old entry, to not duplicate it
    cached_recent(1).checksums(i) = [];
    cached_recent(1).all_vars(i) = [];
    cached_recent(1).badnesses(i,:) = [];
    return;
  end
% else: has nothing in cache, or called to return also the flys structure
end
%clear cached_recent % attempting to write to global before slow SIMION-call, for use while multiple instances in parallel


fly = runsim3(workbench_name, mass,energies, source_point_distribution, variables, counter_local);
% fly = transpose(fly); % turn into row array, to fit the form used in batch runner script_VMI_overview_1 etc.
if isempty(fly) || ~isstruct(fly)
  % Typically when "Warning: No runs selected", i.e. because particles are trapped, going the wrong way or not hitting intended detector.
%   badness = Inf;
  badness = Inf(1, objective_mode); % attempt to handle multiple objectives (DEBUG: improve by whatever logic is needed to tell actual number of objectives)
  return;
end
% if store_result
%   global flys;
%   flys = fly; % DEBUG to make accessible for plotting by call to show_mass_resolution.m from optimize_potentials_loop.m. Not so nice to have enabled when optimizing in parallell though.
% -- not needed, now returned as second argument if nargout>=2
% end

fly_overall_cells = {fly.overall};

% The number of particles missing their detector, initially as array with one column per fly.
% Will be reduced in accordance with objective_mode to have the same size as the badness array.
misses = [fly.misses]; 

more_tolerant_to_misses = ~isempty(strfind(quantity, 'tolmiss'));
    
cost_for_potentials = [];
[m,i_start,i_end] = regexp(quantity, ' pot#(\d+)', 'tokens');
if ~isempty(m)
  % Can use a string like pot#5 to add a cost for potential #5's value.
  %cost_for_potentials = 0.01/1000 * abs(fly(1).potentials(str2double(m{1}))); % An additional 1000V will be as bad as going from 1% to 3% resolution.
  cost_for_potentials = 0.02/1000 * abs(fly(1).potentials(str2double(m{1}))); % An additional 1000V will be as bad as going from 1% to 3% resolution. From 2016-04-11
  %quantity = strrep(quantity, ' pot#',''); % hide to not cause problem with the following processing
  quantity = [quantity(1:i_start-1)  quantity(i_end+1:end)]; % hide to not cause problem with the following processing
end


if isfield(fly(1,1).overall, quantity)
  if ~isempty(strfind(quantity, '_rel'))
    % Note division to convert from relative error (std/avg) into resolving power (avg/std).
    
    %badness = -1./collect_field(fly_overall_cells, quantity, 2); % -Resolving power, better value mapped to more negative badness.
    badness = collect_field(fly_overall_cells, quantity, 2); % Relative error, to use a positive quantity. Giving one array entry per energy
  else
    % Standard deviation, better value is already smaller. Giving one array entry per energy
    badness = collect_field(fly_overall_cells, quantity, 2);
  end
  
elseif ~isempty(strfind(quantity, 'spheres'))
  % Not for optimization, but to show a fine simulated VMI image
  warning('The "spheres" quantity should not be used for optimization.') 
  badness = 1000;
  
else
  % Need computation, not just reading out a precomputed field value.
  
  thetas = fly(1,1).rays(:,strcmp(fly(1,1).ray_columns, 'theta'));
  if fly(1,1).imaging_dimensions == 2
    improved = flipud(linspace(-90, 90, length(thetas))');
    if any(abs(thetas - improved) > 1)
      error('Unexpected angles, rounding error > 1 degree. %s', mat2str(thetas));
    else
      thetas = improved;
    end
  else % 3D
    % TODO improve angle resolution by recomputing
  end
  thetas = thetas * pi/180; % convert to radians
  sin_weights = abs(sin(thetas));

  % This version averaged the energies, without option of objective_mode~=1 
  %sin_weights = repmat(sin_weights' / (sum(sin_weights)*length(fly)), 1, length(fly)); % normalize, to return angle-weighted average per energy, then unweighted average of that across all selected energies
  %rays = collect_field(fly,'rays'); % concatenates the energies, into a (length(fly)*length(thetas))-by-10 matrix
  %y2_rel_sin = sin_weights * (1./rays(:,10)); % |sin theta|-weighted y^2-resolving power (y^2-mean divided by y^2-std)

  sin_weights = sin_weights' / sum(sin_weights); % normalize, to return angle-weighted average per energy, then unweighted average of that across all selected energies
  rays = collect_field(fly,'rays',2); % concatenates the energies, into a length(thetas)-by-(10*length(fly)) matrix
  y2_rel_sin =  1./ ( sin_weights * (1./rays(:,10:10:end)) ); % inverse of: |sin theta|-weighted y^2-resolving power (y^2-mean divided by y^2-std)

  % From the readout here, the list of electron energies doesn't matter
  quantity = strrep(strrep(strrep(strrep(quantity, 'y30_','y_'), 'y20_','y_'), 'y42_','y_'), 'y10_','y_');
  % NOTE WARNING: in case a  new y??_ range is added to optimize_potentials.m but not here in the _objective.m,%
  % it will not be recognized as having electron VMI and will fail to add misses for electron side to badness score, and probaly only optimize for mass (not show second column in Pareto plots)

  % Hide options like ' iso' and ' h3en' from quantity string, to make also the electron-only part of optimize_potentials_objective.m accept it.
  quantity = strrep(strrep(quantity, ' iso',''), ' h3en','');

  switch quantity
    case 'y2_rel sin'
      %badness = -1./y2_rel_sin; % -Resolving power, better value mapped to more negative badness.
      badness = y2_rel_sin; % Relative error, to use a positive quantity, giving one array entry per energy

    case  'y_rel mix pow' % average of 1/'y_rel', 1/'y2_rel', 1/'y2_rel sin', 1/'y2_rel sin' - intended for 2D VMI
      % The divisions convert relative error into resolving power
      y_rel  = collect_field(fly_overall_cells, 'y_rel', 2);
      y2_rel = collect_field(fly_overall_cells, 'y2_rel', 2); % (not sin-weighted like y2_rel_w)

      % Average with double weight on y2_rel_sin, giving one array entry per energy
      %badness = (- 1./y_rel -  1./y2_rel - 2./y2_rel_sin) / 4; % -Resolving power, better value mapped to more negative badness.
      badness = 4 ./ (1./y_rel + 1./y2_rel + 2./y2_rel_sin); % Average as resolving powers, but turn result into relative error to use a positive quantity
    
    case 'y_rel mix' % average of 'y_rel', 'y2_rel', 'y2_rel sin', 'y2_rel sin' - intended for 2D VMI
      % The divisions convert relative error into resolving power
      y_rel  = collect_field(fly_overall_cells, 'y_rel', 2);
      y2_rel = collect_field(fly_overall_cells, 'y2_rel', 2); % (not sin-weighted like y2_rel_w)

      % Average with double weight on y2_rel_sin, giving one array entry per energy
      badness = (y_rel + y2_rel+ 2*y2_rel_sin) / 4; % Average as relative error to use a positive quantity

      if length(badness) < length(energies)
        % This can happen if last particle in a group ends up outside simulation volume without hitting any eletrode.
        % Try lowering Z_GEOMETRY_MAX in Lua script, or enclose the volume by electrodes.
        warning('optimize_potentials_objective:missing_group', 'Missing log summary for %d energies. This happens if last particle in a group ends up outside simulation volume without hitting any eletrode.\n Try lowering Z_GEOMETRY_MAX in Lua script, or enclose the volume by electrodes.', length(energies)-length(badness));
        badness(end+1:length(energies)) = NaN;
      end
      
    otherwise
      % Optimize a combination of relative y-width and mass resolution (requires that both ions and electrons were simulated)
      %e.g. case {'mass', 'y_rel mix + mass', 'y_rel mix + mass3', 'y_rel mix + mass4'}
      if isempty(strfind(quantity, 'mass')) % raise error if 'mass' was not mentioned
        error('Quantity %s not implemented yet...', quantity);
      end
      scaled = ~isempty(strfind(quantity, 'scaled'));
      uncentered = ~isempty(strfind(quantity, 'uncentered'));
      
%       resolutions = show_mass_resolution(fly(1), 0.035, true); % showing figure
      %badness(i) = 1 ./ repmat(resolutions(i).mass_extrapolated_by_spread, 1, length(fly)); % DEBUG currently show_mass_resolution returns one value, after regrouping fly to ignore energy and group by mass.
      badness = NaN(1,length(fly));
%     for i = 1:length(fly)
%       resolutions = show_mass_resolution(fly(i), 3.5, 35, true); %pause; % showing figure
%     end
      if length(badness) < length(energies) % Not sure this check is exhaustive, maybe energies can have a single entry and just different masses. May at laest give warning in typical case.
        % This can happen if last particle in a group ends up outside simulation volume without hitting any eletrode.
        % Try lowering Z_GEOMETRY_MAX in Lua script, or enclose the volume by electrodes.
        warning('optimize_potentials_objective:missing_group', 'Missing log summary for %d energies. This happens if last particle in a group ends up outside simulation volume without hitting any eletrode.\n Try lowering Z_GEOMETRY_MAX in Lua script, or enclose the volume by electrodes.', length(energies)-length(badness));
        badness(end+1:length(energies)) = NaN;
      end

      is_electron = collect_field({fly.source}, 'mass', 2) < 1; % to handle a mixture of ions and electrons flewn
      
      % Since transporting high-energy ions along long drift tube is difficult, and angular distribution not needed, accept some loss without penalty if it is small
      if more_tolerant_to_misses
        % option to accept more misses 2020-02-27, especially when large source like "FWHM5" is used
        to_consider = ~is_electron; % accept some misses on the reflectron side, not in VMI
        has_acceptable_misses = to_consider & misses <= 0.08 * fly(1).source.ion_count;
        misses(has_acceptable_misses) = 0.01 * misses(has_acceptable_misses); % keeping penalty, but making it much smaller if loss is at most 8%
        to_consider = to_consider & ~has_acceptable_misses; % don't apply the following reduction to the same particle(s)
        has_acceptable_misses = to_consider & misses <= 0.15 * fly(1).source.ion_count;
        misses(has_acceptable_misses) = 0.02 * misses(has_acceptable_misses); % keeping penalty, but making it much smaller if loss is at most 15%
        to_consider = to_consider & ~has_acceptable_misses; % don't apply the following reduction to the same particle(s)
        has_acceptable_misses = to_consider & misses <= 0.25 * fly(1).source.ion_count;
        misses(has_acceptable_misses) = 0.2 * misses(has_acceptable_misses); % keeping penalty, but making it smaller if loss is at most 25%
      else % original:
        %misses(~is_electron & misses <= 0.15 * fly(1).source.ion_count) = 0;
        % Changed 2015-05-28 12:38 since outliers and lack of them due to miss instead makes mass resolution readouts uncertain
        has_acceptable_misses = ~is_electron & misses <= 0.15 * fly(1).source.ion_count;
        misses(has_acceptable_misses) = 0.2 * misses(has_acceptable_misses); % keeping penalty, but making it smaller if loss is acceptable (at most 15%).
      end      
      misses_per_fly = misses;

      %[resolutions] = show_mass_resolution(fly(~is_electron), 3.5, 35); % old (until 2015-04-27)
      switch workbench_name
        case {'both8eu', 'both8eu_B2'} % without vertical symmetry, as actually built and wired
          % el.detectorB_x is Lua's detectorB_linecoeff = {k, -1, -d, MCP_centre_y}: k=-1/tan(det.angle=6deg), 
          % d = k*MCP_centre_x - 1*MCP_centre_y; using MCP imported as separate PA (MCP_centre_x=386 as the MCP import Xwb)
          detectorB_x = [-9.5144, -1.0000, 3805.5447, 133.00];
          
        case {'ecyl8e'}
          detectorB_x= [0, 0, -999, -999]; % something invalid, as no ion detector is present in this geomtry

        otherwise
          % TODO: actually as long as the Lua script knows the coefficients, we can now get it from log:
          % detectorB_x = flys(end).geometry.detectorB_linecoeff;
          % flys(end).geometry.detectorB_xSIMION_zMatlab 
          % flys(end).geometry.detector_A_coordinates(3) vertical centre coordinate
          
          warning('optimize_potentials_objective:detectorB_x', 'No detectorB_x value specificly set for workbench in optimize_potentials_objective. Using default.\n May cause particle impacts on detector to be missed.');
          error('No detectorB_x value specificly set for workbench in optimize_potentials_objective. Using default.\n May cause particle impacts on detector to be missed.');
          detectorB_x  = 4.5; %[mm] in between 3.5mm (900mm workbenches) and 5.5mm (1300mm workbenches)
      end
      [resolutions] = show_mass_resolution(fly(~is_electron), detectorB_x , 40); % not showing figure.
      % Note: if multiple masses are given, resolutions contains an entry for each group of two (starting from the heaviest) affect the mass resolutions structure!

      % To not let weighting between mass and y_rel depend on how many kinetic energies included, reshape to one array (mass_basness) for ions mass and another (later) for electrons.
      % This also solves the conversion from resolutions.() where results are per mass (may include multiple energies) rather than per fly.

      if length(resolutions)*2 ~= sum(~is_electron)
        % Didn't get complete result
        mass_badness = NaN(1,sum(~is_electron)/2); % one column per mass-case (each case simulates the pair of masses m & m+1 to determine their separation)
      else
        mass_badness = NaN(1,length(resolutions)); % one column per mass-case (each case simulates the pair of masses m & m+1 to determine their separation)
      end
      for i = 1:length(resolutions)
        if isnan(resolutions(i).mass_extrapolated_by_spacing) % only one mass used
          if length(resolutions(i).mass_used) > 1
            if nargout > 1
              global flys; flys = fly;
              % DEBUG NOTE: raising error here is useful for finding configuration error. But if optimizing over wide range (genetic)
              % it may be more interesting to return a high cost to let optimizer resume with other settings.
              % Nargout=1 when optimizing.
              error('Despite multiple masses used, no .mass_extrapolated_by_spacing result. This suggests the detectorB_x=%s\n in optimize_potentials_objective.m is wrong, or that very bad potentials were given.\n%s', mat2str(detectorB_x), mat2str(fly(1).potentials'));
            end
          end
          if ~isempty(resolutions(i).mass_extrapolated_by_spread) && all(isfinite(resolutions(i).mass_extrapolated_by_spread))
            % use only spread-estimate, and add a penalty (probably smaller than for misses, although having all particles miss may cause this)
            mass_badness(i) = 1 ./ abs(resolutions(i).mass_extrapolated_by_spread) + 1;
          else
            mass_badness(i) = 1000; % no info was available, just put a big penalty
          end
        else % There were two different ion masses
           % Average the two methods. The _by_spacing is a scalar and sensitive to outliers/near-misses that affect the spacing (from average) of the masses,
           % but it is otherwise trusted a bit more. Using average for stabiliity during optimization.
          mass_badness(i) = (mean(1./resolutions(i).mass_extrapolated_by_spread)  +  1/resolutions(i).mass_extrapolated_by_spacing) / 2;
          if abs(diff(resolutions(i).mass_used)) ~= 1
            warning('optimize_potentials_objective:spacing', 'show_mass_resolution called with masses separated by %d u. Suggests many misses! Adding penalty.', diff(resolutions(i).mass_used));
            mass_badness(i) = mass_badness(i) + 10;
          end
        end
        if mass_badness(i) <= 0
          % When masses are not separated, negative mass resolution may be found.
          mass_badness(i) = 500 - 1./mass_badness(i); % Give a large positive number, larger the more negative the mass resolution was
        end
      
        if ~uncentered
          % Add a small penalty in case the image is far from centre of detector
          % The penalty is quadratic with y_mean, up to 40 mm where penalty corresponds to
          % 10 u reduction in resolving power from 200u: (1/190-1/200).
          % However the detector radius is 20 mm, and since the ratio is squared to not care about small deviations, the max penalty was actually 0.25*(1/190-1/200).
          % Changing the division from /40mm to /30mm now to get the squared ratio up to (20/30)^2=0.44 theoretically.
          if ~isnan(resolutions(i).vertical_mean) % 2020-02-26 also vertical axis considered
            mass_badness(i) = mass_badness(i) + (1/190-1/200) * min(1, (resolutions(i).y_mean^2 + resolutions(i).vertical_mean^2)/30^2); % 2020-02-26
            % 2020-02-26 reducing weights for the y_std by dividing not by 40 or 25 but by 50, since a spread is not obviously bad if the positions are still centered and the mass (TOF) resolving power good.
            mass_badness(i) = mass_badness(i) + (1/190-1/200) * min(1, hypot(resolutions(i).y_std, resolutions(i).vertical_std)/50); % 2020-02-26
          else
            %mass_badness(i) = mass_badness(i) + (1/190-1/200) * min(1, resolutions(i).y_mean/40)^2; % 2016
            mass_badness(i) = mass_badness(i) + (1/190-1/200) * min(1, resolutions(i).y_mean/30)^2; % 2020-02-26
            % Also add a corresponding penalty for width of beam, not squared since also small widths are interesting to minimize
            %mass_badness(i) = mass_badness(i) + (1/190-1/200) * min(1, resolutions(i).y_std/40); % 2016
            mass_badness(i) = mass_badness(i) + (1/190-1/200) * min(1, resolutions(i).y_std/50); % 2020-02-26
          end
          
          % on the detector, include vertical_mean in addition to y_mean

          
          % NOTE this cost is quite large now that 7500u is desirable (and >10000 achieved in some geometries). Made option 'uncentered' to allow disabling to find out how detector ought to be moved.
        end
        
        if scaled
          relative_resolving_power = (1/mass_badness(i)) / resolutions(i).mass_used(1);
          % Finally, to improve the averaging of results for multiple masses in different kinetic energy ranges,
          % use a nonlinear function that lets resolving power level-off in case the test mass (desired mass) is resolved.
          % I.e. for 2eV@100u, it should be more important to improve from 80 to 90 u than from 110 u to 120 u.
          %
          % The tanh(( result/aim + 0.1 ).^3) lies in the range 1E-3 to 1.0 and levels off around result/aim=1. result/aim=0.5 ==> 0.21, result/aim=1 ==> 0.87
          % Rescale so that result/aim=1 ==> 1, range is 1.15E.3 to 1.1501.
          %relative_resolving_power = (tanh((relative_resolving_power + 0.1).^3) / tanh(1.1^3));
          relative_resolving_power = (tanh((relative_resolving_power + 0.1)) / tanh(1.1)); % 2015-06-18 14:20 avoid cubing, it made improvements in score beyond relavive power 1.7 completely negligible (not in fifth digit)
          % Alternative without extra score when requirement is exceeded
          % relative_resolving_power = min(1, relative_resolving_power).^3;
          %
          % Convert back to badness (1/resolving power)
          %mass_badness(i) = 1 / (relative_resolving_power * resolutions(i).mass_used(1)); % scaled to a mass
          mass_badness(i) = 1 / (relative_resolving_power * 100); % no longer mass, just arbitrary number (~percentage). This ensures aims at different masses are comparable in the averaging.
        end
      end

      if objective_mode < 3
        % A single output for average mass resolution
        % ( Old, probably never used: if objective_mode == 2, badness(1) = rms(mass_badness); else ... )
        badness = mean(mass_badness); % arithmetic mean of cost function for the masses
        misses = sum(misses_per_fly(~is_electron)); % NOTE: here all non-electrons are considered while only groups of two (from heaviest) affect result of show_mass_resolution()
        objective_index_mass = 1;
      else % objective_mode >= 3
        % New, for Pareto frontier search can return separate objectives even for each mass-case (NOTE: only two cases expected, not tested with one or three...).
        % (For objective_mode == 2 a single mass-objective is still used, presumably using one or several other objective(s) for electron-VMI.)
        badness = mass_badness; % two columns
        objective_index_mass = 1:length(badness);
        misses = misses_per_fly(~is_electron); % NOTE: here all non-electrons are considered while only groups of two (from heaviest) affect result of show_mass_resolution()
        if length(misses) == 2 * objective_index_mass(end)
          % OK, seems to have info for the two masses per test case. Assume that the first are for the higher mass.
          misses = sum(reshape(misses, [2 objective_index_mass(end)])); % reduce to one column per test case, summing misses within each test case
        else
          warning('optpot:unexpected_ion_count', 'Warning: Unexpected number of ion masses used. The number of misses will be averaged among the ion groups.');
          misses = sum(misses)/objective_index_mass(end) * ones(1,objective_index_mass(end)); % distribute the number of misses equally (may be non-integer) between the test cases.
        end
      end
      
      if nargout >= 2
        % When run interactively, print each mass' score
        K = cell2mat({resolutions.kinetic_energy_used}); K = K(1:2:end);
        if scaled
          disp(sprintf('Mass costs (scaled): [%s] =1./[%s %%] for energies %s eV', num2str(mass_badness,' %.6f'), num2str(1./mass_badness,' %.1f'), mat2str(K)))
        else
          disp(sprintf('Mass costs: [%s] =1./[%s u] for energies %s eV', num2str(mass_badness,' %.6f'), num2str(1./mass_badness,' %.1f'), mat2str(K)))
        end
        if nargout >= 3
          % To make accessible e.g. for plotting by call to show_mass_resolution.m from optimize_potentials_loop.m.
          for i = 1:length(fly)
            % Make the detector center coordinate acccessible (for geometries where defined it is the fourth element of detectorB_x)
            % (See also flys(i).source_point_y which is set by runsim3, and for many geometries gives the detectorA's centre coordinate)
            fly(i).detectorB_x = detectorB_x;
          end
        end
      end
  
      if ~isempty(strfind(quantity, 'mass2')) 
        % Double weight for (each) mass resolution -- since 1/50 is OK for VMI and 1/100 is the "scaled" mass score in case the targets are reached without margin.
        badness = badness * 2;
      elseif ~isempty(strfind(quantity, 'mass3'))
        % Triple weight for (each) mass resolution -- since 1/33.333 is decent for VMI (and A-lengh is probably not optimal) but 1/150 expected (1/200u desired) for 1eV ion.
        badness = badness * 3;
      elseif ~isempty(strfind(quantity, 'mass4'))
        % Quadruple weight for (each) mass resolution. NOTE: this is perhaps not reasonable when "scaled" option is used.
        badness = badness * 4;
      end

      
      if ~isempty(strfind(quantity, 'y_rel mix'))
        % Compute also y_rel-mix for electron-VMI, to give additional objectives or an average
        
        y_rel  = collect_field(fly_overall_cells(is_electron), 'y_rel', 2);
        y2_rel = collect_field(fly_overall_cells(is_electron), 'y2_rel', 2); % (not sin-weighted like y2_rel_w)
        %until 2015-12-18 11:59: badness_per_energy = (y_rel + y2_rel+ 2*y2_rel_sin(is_electron)) / 4; % Average quantities with double weight on y2_rel_sin
        %until 2015-12-30 : badness_per_energy = (2*y_rel + y2_rel+ 2*y2_rel_sin(is_electron)) / 5; % Average quantities with double weight on y2_rel_sin and y_rel, to typically avoid high-scoring when y2_rel optimized for lowest energy at the cost of other
        badness_per_energy = (y_rel + y2_rel+ y2_rel_sin(is_electron)) / 3; % Average quantities with equal weight. In recent improvement of Lua code, it was discovered that no y-center coordinate was used (thus falsly good resolution when detector off-centre by 87 mm). Now it seems OK to average with equal weight
        
        % TODO IMPROVEMENT in y2_rel_sin(~is_electron) VMI of ions can also be favoured (perhaps with some small coefficient)
        
        if objective_index_mass(end) < 4
          % Append a single objective representing all electron energies
          objective_index_electron = objective_index_mass(end)+1;
          misses(objective_index_electron) = sum(misses_per_fly(is_electron));

          %until 2015-12-18 11:59: badness(2) = mean(badness_per_energy); % Arithmetic average among energies
          %until 2016-01-06 15:50: badness(2) = sqrt(mean(badness_per_energy.^2)); % 2-norm average among energies, to put more weight on worst energy, striving to get all OK rather than one very good and others worse
          badness(objective_index_electron) = mean(badness_per_energy.^3)^(1/3); % 3-norm average among energies, to put even more weight on worst energy, striving to get all OK rather than one very good and others worse

          
        else % if objective_index_mass >= 4
          % Separate objectives for each electron energy, possibly useful for Pareto front search.
          objective_index_electron = objective_index_mass(end) + 1:length(badness_per_energy);
          misses(objective_index_electron) = misses_per_fly(is_electron); % (these array sizes should match)
          badness(objective_index_electron) = badness_per_energy;     
          
        end
        
        if nargout >= 2
          % When run interactively, print e-VMI score
          K_electrons = collect_field({fly(is_electron).source}, 'energy', 2);
          disp(sprintf('Electron-VMI costs: [%s] =1./[%s] for energies %s eV', num2str(badness_per_energy,' %.6f'), num2str(1./badness_per_energy,' %.1f'), mat2str(K_electrons)))
        end
      else
        % Ignore electrons, just consider mass resolution
        % ignore electron misses
      end

%     otherwise
%       error('Quantity %s not implemented yet...', quantity);
  end
  
end

which_NaN = isnan(badness);
if any(which_NaN)
  % To avoid NaN or Inf, replace any NaN result with a normal result plus the penalty for misses.
  misses(which_NaN) = max(misses(which_NaN), 1); % if NaN without miss, treat as 1 miss anyway
  if any(~which_NaN)
    badness(which_NaN) = max([NaN, badness(~which_NaN)]); % NaN remains if there was no miss-free energy
  end
end

if any(misses > 0)
  % Add a large number to mark cases where particles miss detector as bad,
  % but still allow optimizer to tell the difference between cases with few or many misses.
  which_missed = find(misses > 0);
  if more_tolerant_to_misses
    penalty = 20 * ones(1,length(which_missed));
  else
    penalty = 100 * ones(1,length(which_missed));
  end
  % Added 2015-05-28 12:38 to reduce the fixed-term penalty when few outliers, regardless of on electron or ion side. 
  % (Resolution results will not be reliable but optimizer will get a strong nudge to approach zero-miss configuration.)
  % If not more than 5% of the particles of any type missed, reduce the big penalty term
  % to help optimizers realize that this is near a good configuration
  has_acceptable_misses = misses(which_missed) <= 0.05 * fly(which_missed(1)).source.ion_count;
  penalty(has_acceptable_misses) = 1; % instead of 100 (the fixed penalty term, the actual number of misses is used on the next line to compute basness)
  
  if more_tolerant_to_misses % option 2020-02-27
    has_acceptable_misses = misses(which_missed) <= 0.20 * fly(which_missed(1)).source.ion_count; % at most 20%
    badness(which_missed(has_acceptable_misses)) = badness(which_missed(has_acceptable_misses)) + 0.1*sum(max(abs(badness(isfinite(badness))))) + penalty(has_acceptable_misses) + misses(which_missed(has_acceptable_misses));
    badness(which_missed(~has_acceptable_misses)) = badness(which_missed(~has_acceptable_misses)) + sum(max(abs(badness(isfinite(badness))))) + penalty(~has_acceptable_misses) + 2*misses(which_missed(~has_acceptable_misses));
  else % original
    badness(which_missed) = badness(which_missed) + sum(max(abs(badness(isfinite(badness))))) + penalty + 2*misses(which_missed); % (before 2016-04-05 was factor 10* before misses(which_missed))
  end
  
  if objective_mode > 1
    % Objectives won't be averaged. To reduce the number of instances where the Pareto optimizer
    % selects settings that miss e.g. the highest energy completely just to favour the lower energies,
    % add a penalty also to other energies than the one that missed. From 2016-04-06.
    if all(has_acceptable_misses) % 2020-02-27
      badness(~which_missed) = badness(~which_missed) + 0.025;
    else % original
      badness(~which_missed) = badness(~which_missed) + 0.1;
    end
  end
  
end
% 2016-01-12 Maybe the extremely high badness when seveal misses are not so helpful (e.g. in Pareto plot)?
% Reduce the growth rate for badness beyond a threshold (per objective). As the fitness sorting is actally by rank in normal genetic, it doesn't matter as long as the mapping is monotonous.
threshold = 5;
badness(badness > threshold) = threshold + 0.01 * (badness(badness > threshold)-threshold);

% TODO FIXME: the costs_before_averaging output may need to be produced earlier, since if non-Pareto the average is taken earlier no than in old version.
% This may cause poblem when running to show plots and get output info as second and third argument.
costs_before_averaging = badness;

% If multiple groups are flewn (e.g. different energies), how to combine their results? 
% switch objective_mode
%   case 1 % Arithmetic average. Probably always used until 2016-01-11
%     badness = mean(badness);
%   % Old intent but normally not used
%   %case 0 % Don't combine, output a vector (which optimizer will optimize 2-norm of).
%   %  % This gives the optimizer more info than when a scalar is returned.
%   %case 2 % 2-norm (like 0 but performed here)
%   %  badness = rms(badness);
% 
%   case 2
%     % For Pareto frontier search, return one column per objective.
%       ... This is already handled above when populating the badness and misses arrays.
% 
%   otherwise
%     error('Unsupported objective_mode method %d', objective_mode);
% end


% if ~isempty(strfind(quantity, 'pot')) && strcmp(variables{end-1},'potentials')
%   
%   TODO cost based on max |potential| in  variables{end}, and maybe some for e-drift tube
%   
% end

badness(~isfinite(badness)) = Inf; % Change a possibly remaining NaN into +Inf to not stop some optimizers
%badness(~isfinite(badness)) = 1E6; % Change a possibly remaining NaN/Inf into +large number to not stop any optimizer (but gives bogus printouts of really bad "optima" when cause most likely was error in workbench-Lua configuration)

if ~isempty(cost_for_potentials)
  % Add cost_for_potentials to each index. (Having it as separate objective didn't seem to improve Pareto results in the intended direction, and for non-Pareto it doesn't matter how it enters average.)
  badness = badness + cost_for_potentials;
end

% DEBUG show potentials and fitness values (before averaging in case objective_mode==1)
index = find(strcmp(variables,'potentials')) + 1;
if ~isempty(index)
%   disp(sprintf('%s\t: %.6g', num2str(variables{index}(:)', '%8.2f  '), badness))
%  disp(sprintf('%s : %.6g', sprintf('%8.2f  ', variables{index}(:)'), badness)) % with 0.01 V resolution
  if length(badness) == 1
    disp(sprintf('%s : %.6g', sprintf('%5.0f ', variables{index}(:)'), badness)) % with 1 V resolution
  else % for Pareto search using multi-objective function
    disp(sprintf('%s : %s', sprintf('%5.0f ', variables{index}(:)'), sprintf('%.5f ', badness)))
  end
end

if objective_mode <= 1 % < 2
  % Return a single output cost, averaging the mass- and electron- cost functions
  % This is an arithmetic mean.
  badness = mean(badness); 
else
  
end

if nargout >= 3
  % To make accessible e.g. for plotting by call to show_mass_resolution.m from optimize_potentials_loop.m.
  flys = fly; 
  %   size(fly) % DEBUG
end


if nargout == 1 
  %global cached_recent % attempting to read from global a second time after slow SIMION-call, for use while multiple instances in parallel
  if isempty(cached_recent)
    % Initialize cache
    %disp(sprintf('Thread %.8f, counter %d. First cached %.16g with %.8f.', optpot_thread_id, counter_local, checksum, badness)); % DEBUG
    cached_recent(1) = struct('checksums', {checksum}, 'all_vars',{{all_vars}}, 'badnesses',{badness});
  else
    % Append to cache
    %disp(sprintf('Thread %.8f, counter %d (%d). Appending %.16g as #%d with %.8f.', optpot_thread_id, counter_local, counter_local_raw, checksum, length(cached_recent(1).checksums)+1, badness)); % DEBUG
    cached_recent(1).checksums(end+1,1) = checksum;
    cached_recent(1).all_vars{end+1,1} = all_vars;
    cached_recent(1).badnesses(end+1,:) = badness;

    % if size(cached_recent(1).checksums,1) > 30 % unless more than 15 individuals in genetic, this is is sufficient even if assigned to wrong thread in every second generation
    if size(cached_recent(1).checksums,1) > 60 % unless more than 30 individuals in genetic, this is is sufficient even if assigned to wrong thread in every second generation
      % Limit length of cache buffer by removing oldest entry
      %disp(sprintf('Thread %.8f. Dropping oldest %.16g.', optpot_thread_id, cached_recent(1).checksums(1))); % DEBUG
      cached_recent(1).checksums(1) = [];
      cached_recent(1).all_vars(1) = [];
      cached_recent(1).badnesses(1,:) = [];
    end

  end
else
  % Auto-clear cache when called to show resuls (nargout>1)
  cached_recent=[];
end
