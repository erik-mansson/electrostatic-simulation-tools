% Assess mass resolution and/or related aspects 
% for one or two given flys-structures (from readsim3).
%
% The highest mass, and the second-to-highest mass will be used. Lower masses ignored.
% If only simulation with one mass are given, some output will be NaN.
%
% PARAMETERS
%  flys          Array of fly-structs from readsim3
%  z_detector    [mm] ignore partices impacting at z-coordinates which are not near z_detector. (Default: 4 mm)
%                This can either be a scalar, z_coordinate, or it can be a four-element array specifying
%                [k1 k2 -d, y_center] for the line equation  k1*z + k2*y - d = 0 (scalar means k1=1,k2=0,d=z_detector).
%  max_r         [mm] (Default: Inf)
%                If z_detector is scalar: ignore partices impacting the detector at larger radius than this. 
%                But if z_detector is an array specifying a line equation, the max_r check is not peformed,
%                setting a finite max_r then just fixes the abscissa range of some diagrams.
%                
%  show_figure   Show figure and print to command window? (Default: false if output arguments received, true otherwise)
%  annotation_suffix  String to append to central annotation
%
% RETURNS
%  resolutions   Struct with the fields
%                mass_extrapolated_by_spacing:   separation_rel_spread * the_lower_mass [u], or NaN if only one mass.
%                separation_rel_spread:          (TOF(mass2) - TOF(mass1)) / mean(t_spread), or NaN if only one mass.
%                mass_extrapolated_by_spread:     TOF/2/spread which is an estimate of max mass wher 1u resolved (1-by-2 or scalar).
%                t_spread:                       the spread of TOF for each mass (1-by-2 or scalar).
%                y_mean:                         average of y-coordinates (scalar)
%                y_std:                          standard deviation of y-coordinates (scalar)
%                type_of_spread:                 Currently '5% to 95% IQR', meaning interquantile range from 5% to 95%.
%                mass_used:                      The mass for which the simulation was made [u] (1-by-2 or scalar).
%                kinetic_energy_used:            The kinetic energy for which the simulation was made [u] (1-by-2 or scalar).
%                max_r_used:                     [mm] the max_r input argument.
% geometric_coefficients
%                Parameters describing how x0, y0, and y correlate to the TOF-deviation.
%                One column per input flys.
%                For coeff_y, the rows are the coefficients for constant, y^2 and y^4 of an even-order-polynomial fit.
%                For coeff_y0, the rows are the coefficients for constant, y0^2 and y0^4 of an even-order-polynomial fit.
%                Q_r0t is a fit of the dimensionless first derivative v1 * dt/d(r0)
%                Q_z is a fit of the dimensionless first derivative v1 * dt/d(z_0)
%                where v1 is the speed at extractor due to extraction work (from nominal source point).
%                      v1 = sqrt(2 * nominal_extraction_work / mass);
%
% i2             The indices in flys and geometric_coefficients that belong to the higher mass.
% i1             The indices in flys and geometric_coefficients that belong to the lower mass.
% handles        Struct with some handles to parts of figure
%
% SEE ALSO readsim3 runsim3
function [resolutions, geometric_coefficients, i2, i1, handles] = show_mass_resolution(flys, z_detector, max_r, show_figure, annotation_suffix)
constants;
handles = struct();

show_x_instead_of_r = false; % default
show_x_instead_of_r = true; % OPTION (will use other filename)


% TODO: consider changing output, to return one mass_extrapolated_by_spread-value per flys-entry, which could be different kinetic energies for the same mass. But can be emulated by calling with only flys for desired energy, e.g. flys(1).

if nargin < 2
  % Originally
  %z_detector = 4;%3.5; %[mm] assume detector is on B-side, near edge of workbench, without rotation
  z_margin = 4;
  error('Detector B position or coefficients not defined. TODO could use from Lua log file now. detectorB_linecoeff, detectorB_xSIMION_zMatlab, detector_A_coordinates');
  % TODO: actually as long as the Lua script knows the coefficients, we can now get it from log:
  % detectorB_x = flys(end).geometry.detectorB_linecoeff;
  % flys(end).geometry.detectorB_xSIMION_zMatlab 
  % flys(end).geometry.detector_A_coordinates(3) vertical centre coordinate
else
  z_margin = 1.5; %[mm], assume more accurate if z_detector parameter was given
end
if nargin < 3
  max_r = Inf;
end

if nargin < 4
  show_figure = nargout == 0;
end
flys = transpose(flys(:)); % turn into row vector
if nargin < 5
  annotation_suffix = '';
end

if flys(1).source.mass < 1 && flys(end).source.mass < 1 
  error('Trying to check mass resolution using simulation with only electrons.');
end

if length(flys(1).source.spacing) == 3 % 3D Gaussian
  source_FWHM_y_str = sprintf('W%.1f', flys(1).source.spacing(2) * (2*sqrt(2*log(2)))); %[mm]
  annotation_suffix = [annotation_suffix sprintf('. Source FWHM: flight axis %.1f, optical axis y %.1f, vertical %.1f [mm].', ...
            flys(1).source.spacing(1) * (2*sqrt(2*log(2))), flys(1).source.spacing(2) * (2*sqrt(2*log(2))), flys(1).source.spacing(3) * (2*sqrt(2*log(2))) )]; %[mm]
else
  source_FWHM_y_str = '';
end

% Get general coordinate info
all_trajectories = {flys.trajectories};
if isfinite(max_r) || ~isempty(z_detector)
%   all_which = abs(collect_field(all_trajectories, 'y')) <= max_r;
  if length(z_detector) == 1
    % For cylinder symmetry about z-axis
    all_which = (abs(collect_field(all_trajectories, 'y')) <= max_r) & (abs(collect_field(all_trajectories, 'z')-z_detector) <= z_margin);
  else
    % Detector surface is a line when projected onto (z,y)-plane. Check that deviation from that line is smaller than z_margin
    all_which = abs( z_detector(1)*collect_field(all_trajectories, 'z') ...
                   + z_detector(2)*collect_field(all_trajectories, 'y') + z_detector(3) ) <= z_margin;
    % ignoring max_r
  end
else
  all_which = ':';
end
all_z0 = collect_field(all_trajectories, 'z0');
all_z0 = all_z0(all_which);
if isfield(flys(1).source,'source_point_z')
  nominal_z0 = flys(1).source.source_point_z;
else
  nominal_z0 = mean(all_z0); % approximate the nomininal source point using all given flys, since readsim3 doesn't provide a value for it
end
if isfield(flys(1).source,'source_point_y')
  nominal_y0 = flys(1).source.source_point_y;
else
  nominal_y0 = 0; % assume cylinder symmetric nominal case
  disp('No source_point_y given (e.g. when read directly by readsim3.m instead of optimize_potentials.m).');
end

if length(z_detector) == 1 || length(z_detector) == 3
  % Detector is centered on z-axis (or fourth parameter not given, so coordinate remapping not possible)
  nominal_y = 0;
  nominal_z = NaN;
  remap = @(x,y) y;
  abscissa_y  = (0:0.025:1) * ceil( max(abs(collect_field(all_trajectories, 'y'))) / 0.2)*0.2; %[mm] (currently ignoring max_r filter here, will make it visual what range was skipped (no points there))
  abscissa_y0 = (0:0.025:1) * ceil( max(abs(collect_field(all_trajectories, 'y0')-nominal_y0)) / 0.5)*0.5; %[mm] (currently ignoring max_r filter here, should be OK)
else
  % Determine detector position (not perpendicular to z(SIMION x) axis in reflectron geometries
  nominal_y = z_detector(4);
  if z_detector(1) ~= 0 % normal case, the angle is nonzero (detector surface normal is not parallel to the z-axis)
    nominal_z = (-z_detector(2)*nominal_y - z_detector(3)) / z_detector(1);
  else %special case when detector normal is along z axis (in this case the roles of y and z axes become somewhat reversed to what is assumed mostly)
    nominal_z = z_detector(4);
    nominal_y = (-z_detector(1)*nominal_z - z_detector(3)) / z_detector(2);
  end
  % Project onto a vector along the detector surface (it is rotated 90 degrees wrt. the normal vector in z_detector(1:2))
  %remap = @(z,y) (z_detector(2)*(z-nominal_z) - z_detector(1)*(y-nominal_y)) / sqrt(z_detector(2)^2 + z_detector(1)^2);
  remap = @(z,y) (z_detector(2)*(z-nominal_z) - z_detector(1)*(y-nominal_y)) / norm(z_detector(1:2));
  % NOTE: this 1D projection still ignores the x (third) coordinate, not normally used in simulations. Could use proper [dx;dy;0] - [z_detector(1);z_detector(2);0]*[z_detector(1) z_detector(2) 0]*[dz;dy;dx]/norm projection to get the remaining vector in detector plane

  % Show negative y and y0 too, since symmetry cannot be assumed
  if ~isempty(all_trajectories{1}.t)
    abscissa_y  = (-1:0.025:1) * ceil( max(abs(remap(collect_field(all_trajectories, 'z'),collect_field(all_trajectories, 'y')))) / 0.2)*0.2; %[mm] (currently ignoring max_r filter here, will make it visual what range was skipped (no points there))
    abscissa_y0 = (-1:0.025:1) * ceil( max(abs(collect_field(all_trajectories, 'y0')-nominal_y0)) / 0.5)*0.5; %[mm] (currently ignoring max_r filter here, should be OK)
  else % No hits at all
    abscissa_y = (-1:0.025:1) * 25;
    abscissa_y0 = (-1:0.025:1) * 25;
  end
end

full_dz0_range = (minmax(all_z0') - nominal_z0); %[ns]
if diff(full_dz0_range) == 0
  full_dz0_range = full_dz0_range + [-1 1]; %[ns] ensure some nonszero
end
clear all_z0 all_which;

% Added 2015-05-28 20:44, the sign shown for Q_z was wrong (since detector on negative z-side used)!
% When the drift tube is too long with respect to extraction region and b & f parameters, Q_z should be positive.
if nominal_z0 > z_detector
  z_sign = -1; % using detector on negative z-side (if z=0 in nominal source point)
else
  z_sign = 1;
end

% Mass for each fly (particle type)
masses = collect_field({flys.source}, 'mass', 2);

% These variables have one entry per fly, entries for non-selected particles (e.g. electrons) will be left containing NaN
Q_z = NaN(1,length(masses));
Q_r0t = Q_z;
coeff_y0 = NaN(3,length(masses)); % [constant offset; y0^2-coefficient; y0^4-coefficient], one column per fly
coeff_y  = coeff_y0; % [constant offset; y^2-coefficient; y^4-coefficient], one column per fly
coeff_r  = coeff_y0; % [constant offset;r^2-coefficient; r^4-coefficient], one column per fly
coeff1_x  = NaN(2,length(masses)); % [constant offset; x-coefficient], one column per fly
coeff1_y  = NaN(2,length(masses)); % [constant offset; y-coefficient], one column per fly
ratio_wx_wy = NaN(1,length(masses)); % std(x)/std(y)
h_lines = Q_z;

% Initialize, since in rare cases no output was assigned otherwise (probably when no particle hit the detector)
resolutions(1,1) = struct('mass_extrapolated_by_spacing', NaN, 'separation_rel_spread', NaN, 'mass_extrapolated_by_spread', NaN, ...
  't_spread', NaN, 'y_mean', NaN, 'y_std', NaN, 'vertical_mean', NaN, 'vertical_std', NaN, ...
  'type_of_spread', '', 'mass_used', [], 'kinetic_energy_used', [], 'max_r_used', max_r, 'z_det_used', z_detector);

colors = 'kbgcrmy'; % for different flys of same mass (e.g. if more than one kinetic energy)
for ri = 1:(sum(masses(:)>1)/2) % result index, counts the pairs of masses (the entries in the results struct array)
  % Get a pair of masses to analyse resolution for

  mass2 = max(masses); % the heaviest mass
  i2 = masses == mass2; % which of the flys has the heviest mass
  mass1 = max(masses(~i2)); % the second-to-heaviest mass, if any
  if isempty(mass1)
    i1 = [];
  else
    i1 = find(masses == mass1);
  end
  i2 = find(i2); % NOTE: this can be more than one index, of several flys for each mass.
  if length(i2) > 1 && flys(i2(1)).source.energy ~= flys(i2(end)).source.energy
    warning('Warning: multiple flys with different energies for mass %.0f. Some output may be incorrect.', mass2); % IMPROVEMENT: decide if allowed or not
  end
  if length(i1) > 1 && flys(i1(1)).source.energy ~= flys(i1(end)).source.energy
    warning('Warning: multiple flys with different energies for mass %.0f. Some output may be incorrect.', mass1);
  end

  if min([mass2 mass1]) < 1
    if ri > 1
      % Already got the first mass pair. Stop if only electrons remain,
      % and don't need to try to adding an entry for single mass in case an odd number of masses present. 
      break; 
    else
      error('show_mass_resolution:electron', 'show_mass_resolution called with particles of mass < 1 u.');
    end
  elseif mass2 - mass1 ~= 1
    warning('show_mass_resolution:spacing', 'show_mass_resolution called with masses separated by %d u. The result should not be taken as max mass where 1 u is resolved!', mass2-mass1);
  end
  % TODO give error also if charge differs between the two masses

  % Pick the lens info for the side the particle would have been heading towards
  if flys(i2(1)).lens_parameters_B.U_ext < 0
    lens = flys(i2(1)).lens_parameters_A;
  else
    lens = flys(i2(1)).lens_parameters_B;
  end


  %% Quantify TOF-spread due to source point, for each given flys
  if show_figure
    figure(8); clf;
    %set(8,'Position', [1 31 1600 794]); % corresponds to fullscreen on laptop monitor
    set(8,'Position', [4 23 1276 794]); % on office monitor, the width gets truncated to this (which can be used always to get a consistent size on all monitors)
  end
  counter = 0;
  previous_mass = -Inf;

  % One resolutions-entry per set [i2 and possibly i1] of masses
  resolutions(ri,1) = struct('mass_extrapolated_by_spacing', NaN, 'separation_rel_spread', NaN, 'mass_extrapolated_by_spread', NaN, ...
    't_spread', NaN, 'y_mean', NaN, 'y_std', NaN, 'vertical_mean', NaN, 'vertical_std', NaN, ...
    'type_of_spread', '', 'mass_used', mass2, 'kinetic_energy_used', flys(i2(1)).source.energy, 'max_r_used', max_r, 'z_det_used', z_detector);

  for i = [i2 i1]
    if flys(i).source.mass == mass2
      style = 'o'; % circles for the higher mass
      row = 1;
    else
      style = '.'; % dots for the lower mass, if any
      row = 2;
    end
    style = [style colors(mod(counter, length(colors))+1)];
    counter = counter + 1;

    % v1 = speed at extractor due to extraction work (from nominal source point) = sqrt(2 * nominal_extraction_work / mass)
    v1 = sqrt(2 * abs(flys(i).source.charge*eV * lens.U_ext) / (flys(i).source.mass*atomic_mass_unit)); % [m/s]

    if flys(i).source.mass ~= previous_mass
      % Once for each mass
      trajectories_for_mass = {flys(i).trajectories}; % for this mass
      t_for_mass = collect_field(trajectories_for_mass,'t');
      %flys.source.z0 TODO
      z0_for_mass = collect_field(trajectories_for_mass,'z0');
      if isfinite(max_r) || ~isempty(z_detector)
        if length(z_detector) == 1
          % For cylinder symmetry about z-axis
          which_for_mass = abs(collect_field(trajectories_for_mass,'z') - z_detector) <= z_margin;

          %t_axis_range = minmax(t_for_mass(which_for_mass)'); % before filtering by |y| <= max_r
          which_for_mass = which_for_mass & (abs(collect_field(trajectories_for_mass,'y')) <= max_r);
        else
          % Detector surface is a line when projected onto (z,y)-plane
          which_for_mass = abs( z_detector(1)*collect_field(trajectories_for_mass, 'z') ...
                              + z_detector(2)*collect_field(trajectories_for_mass, 'y') + z_detector(3) ) < z_margin;
        end
        t_for_mass  = t_for_mass(which_for_mass);
        z0_for_mass = z0_for_mass(which_for_mass);
      else
        which_for_mass = ':';
      end
      if isempty(t_for_mass)
        warning('show_mass_resolution:nothing', 'No particle hits accepted for mass %d u.', flys(i).source.mass);
        % resolutions cannot be determined!
        %return;
        break;
      end
        
      t_axis_range = minmax(t_for_mass'); % after filtering
      if diff(t_axis_range) == 0
        t_axis_range = t_axis_range + [-1 1]; %[ns] ensure the axis has a nonzero range
      end
      [intensity,bins] = hist(t_for_mass, 40);
      t_std = std(t_for_mass); %[ns]
      z0_std = std(z0_for_mass);

      % To require rather clean peaks for what we call mass resolution, use 5%--95% interquantile range
      % rather than standard deviation or FWHM (FWHM also difficult on noisy histogram).
      % 5%--95% seems visually more relevant than 10%--90%.
      t_IQR = diff(quantile(t_for_mass, [0.05 0.95]));
      resolutions(ri).type_of_spread = '5% to 95% IQR';

      % t_IQR / t_std is about 3.28 on one example file (for 5-95% range), t_IQR/fwhmi(bins, intensity) is about 1.4 on that example.

      if row == 1 %mass2
        resolutions(ri).t_spread = t_IQR;
      elseif row == 2 % mass1 (< mass2)
        % Now, t_mean_for_mass and t_spread still refer to the previously used mass: mass2, while mean(t_for_mass) and t_IQR refer to mass1
        resolutions(ri).t_spread = [resolutions(ri).t_spread t_IQR];  % prepend, so both the mass1- and mass2-based values are returned if flys for more than one mass was given
        resolutions(ri).separation_rel_spread = (t_mean_for_mass - mean(t_for_mass)) / mean(resolutions(ri).t_spread);% use mean of spread of the two masses as result
        if mass2 - mass1 ~= 1
          % When warning 'show_mass_resolution:spacing'. To ensure automatic optimization does not get good score, divide resolution by mass separation (1 in the normal case)
          resolutions(ri).separation_rel_spread = resolutions(ri).separation_rel_spread / max([1 abs(mass2-mass1)]);
        end

        % TODO: figure out if the relative factor at current mass should be squared or square-rooted
        % to better estimate at what other mass the separation_rel_spread would reach zero. May depend on whether
        % kinetic energy is assumed constant or decreasing for higher mass...
        resolutions(ri).mass_extrapolated_by_spacing = mass1 * resolutions(ri).separation_rel_spread;
        resolutions(ri).mass_used = [mass1 mass2];
        resolutions(ri).kinetic_energy_used = [flys(i1).source.energy flys(i2).source.energy];

        % TODO: check that lower quantile of high mass > higher quantile for low mass. This should be better when outliers (neglecting <5% outliers but not if more).

      end

      t_mean_for_mass = mean(t_for_mass);
      % M/deltaM is approximately T/(2*deltaT). Here shown using delta=5-to-95% percentile range instead of FWHM or std
      mass_resolution_by_IQR = t_mean_for_mass/2/t_IQR;

      if isnan(resolutions(ri).mass_extrapolated_by_spread) % mass2
        resolutions(ri).mass_extrapolated_by_spread = mass_resolution_by_IQR;
      else % mass1 (< mass2)
        resolutions(ri).mass_extrapolated_by_spread = [mass_resolution_by_IQR resolutions(ri).mass_extrapolated_by_spread]; % prepend, so both the mass1- and mass2-based values are returned
      end


      if show_figure
        subplot(2,5,(row-1)*5 + 1);
        barh(bins/1E3,intensity);
        title(sprintf('%.0f u: 5%%IQR=%.0f ns, \\sigma=%.0f ns\n (|Q_z|~%s), t/2IQR=%.0f', flys(i).source.mass, ...
          t_IQR, t_std, texformat_g10(v1 * t_std*1E-9/(z0_std*1E-3)), mass_resolution_by_IQR));
        ylabel('t / \mus')
        axis tight
  %       t_axis_range = ylim();
        if ~isempty(t_axis_range) && all(isfinite(t_axis_range)) && diff(t_axis_range)>0
          ylim(t_axis_range/1E3);
        end
      end
    end

    previous_mass = flys(i).source.mass;

    traj = flys(i).trajectories;
    
    if isfield(flys(end),'geometry') && ~isempty(flys(end).geometry)
      % Not sure it is needed, but for both8eu which is not centered vertically, subtract the offset so vertical coordinate is centered as it was in older models
      traj.x = traj.x - flys(end).geometry.detector_A_coordinates(end);
      traj.x0 = traj.x0 - flys(end).geometry.detector_A_coordinates(end);
    end
    
    if isfinite(max_r) || ~isempty(z_detector)
      if length(z_detector) == 1
        % For cylinder symmetry about z-axis
        which = abs(traj.y) <= max_r & abs(traj.z - z_detector) <= z_margin;
      else
        % Detector surface is a line when projected onto (z,y)-plane
        which = abs( z_detector(1)*traj.z + z_detector(2)*traj.y + z_detector(3) ) < z_margin;
      end
    else
      which = ':';
    end
    t = traj.t(which);
    %fprintf('Got %d of %d ions of mass %.1f u, energy %.5f eV \n', length(t), length(traj.mass), traj.mass(1), traj.energies(1)) % DEBUG
    if length(z_detector) < 4
      y = traj.y(which);
    else
      % Project onto detector basis vector (along detector surface)
      y = remap(traj.z(which), traj.y(which));
    end

    % Prepare to let optimize_potentials_objective.m add a small penalty if the ions are not centered on the detector
    if isnan(resolutions(ri).y_mean)
      resolutions(ri).y_mean = mean(y);
      resolutions(ri).y_std = std(y);
    else
      resolutions(ri).y_mean = mean([resolutions(ri).y_mean mean(y)]); % average the two masses (assumes this won't be run for more than two masses)
      resolutions(ri).y_std = mean([resolutions(ri).y_std std(y)]); % average the two masses (assumes this won't be run for more than two masses)
    end

  %  scaled_dt = v1 * (traj.t-mean(traj.t))*1E-9; %[m] % use mean TOF per fly, i.e. as if kinetic energy-dependent TOF-center-offset could be ignored
    scaled_dt = v1 * (t - t_mean_for_mass)*1E-9; %[m] % alerternative, using same mean for all flys of same mass, more relevant for mass resolution where structure within peak is not accessible.

    % TODO check what we plotted in Laksman spectrometer paper
    
    if length(t) > 4 && length(y) > 4
      coeff_y(:,i) = ( (repmat(y,1,3) .^ repmat([0 2 4],length(y),1)) \ (t - t_mean_for_mass) )';
      coeff1_y(:,i) = ( (repmat(y,1,2) .^ repmat([0 1],length(y),1)) \ (t - t_mean_for_mass) )';
    else % too few hits to fit three coefficients
      coeff_y(:,i) = NaN;
      coeff1_y(:,i) = NaN;
    end
    if show_figure && ~isempty(t)
      subplot(2,5,(row-1)*5 + 2);
      if length(z_detector) == 1
        plot(abs(y), t - t_mean_for_mass, style);
      else
        plot(y, t - t_mean_for_mass, style);
      end
      hold on;
      h_lines(i) = plot(abscissa_y, polyval([coeff_y(end,i) 0 coeff_y(2,i) 0 coeff_y(1,i)], abscissa_y),  ['-' style(end)]);
      plot(abscissa_y, polyval(coeff1_y([2 1],i), abscissa_y),  ['--' style(end)]);
      ylabel('t - <t> / ns')
      if ~isempty(t_axis_range) && all(isfinite(t_axis_range)) && diff(t_axis_range)>0
        ylim(t_axis_range - t_mean_for_mass);
      end
      if resolutions(ri,1).kinetic_energy_used <= .2 %[eV]
        standard_det_halfrange = 5; %[mm]
      elseif resolutions(ri,1).kinetic_energy_used <= 1 %[eV]
        standard_det_halfrange = 10; %[mm]
      else
        %standard_det_halfrange = 30; %[mm]
        standard_det_halfrange = 22; %[mm]  detector is only 40 mm in diameter, thus show 44 rather than 60 mm diameter
      end
      if length(z_detector) == 1 % (no longer used for interesting geometries like both8e and both8eu)
        xlabel('refl.horiz. |y| / mm')
        if diff(minmax(y(:)')) < 2 * standard_det_halfrange
          xlim( mean(minmax(y(:)')) + [-1 1]*standard_det_halfrange )
%         xlim([0 max([ceil(abscissa_y(end)) 30])]) % use fixed axis range unless very broad spread of data
        else
          if diff(minmax(y(:)')) < 6 * standard_det_halfrange
            xlim( mean(minmax(y(:)')) + [-3 3]*standard_det_halfrange )
          else % (old) show full range
            xlim([0 abscissa_y(end)])
          end
        end
      else % used:
        xlabel('refl.horiz. y / mm')
        %if diff(minmax(y(:)')) < 2 * standard_det_halfrange
        %  xlim( mean(minmax(y(:)')) + [-1 1]*standard_det_halfrange )
        % %xlim([-1 1] .* max([ceil(-abscissa_y(1)) ceil(abscissa_y(end))], [30 30])) % use fixed axis range unless very broad spread of data
        if max(abs(y(:))) <= 1 * standard_det_halfrange && ~isfinite(max_r)
           xlim( 0 + [-1 1]*standard_det_halfrange );
        else
          %if diff(minmax(y(:)')) < 6 * standard_det_halfrange
            %xlim( mean(minmax(y(:)')) + [-3 3]*standard_det_halfrange )
          %else % (old) show full range 
          %  xlim([abscissa_y(1) abscissa_y(end)])
          %end
          %xlim( mean(minmax(y(:)')) + [-22 22] ) % as the high-energy range (44 mm diameter is a bit more than the detector)
          if isfinite(max_r)
            xlim([-max_r max_r]);
          else
            xlim([-22 22]); % as the high-energy range (44 mm diameter is a bit more than the detector)
          end
        end
      end

      title(sprintf('Misses: %d/%d here (%d for all %d particles)', flys(i).misses, flys(i).source.ion_count, sum(collect_field(flys,'misses')), length(flys) ));
    end
    
    % BEGIN x
    x = traj.x(which); % this is the vertical axis
%     x = traj.x0(which); % OPTION: this is the vertical axis at source. Try again if simulated with larger spread
%     x = x + (t-t_mean_for_mass) % DEBUG
    r = hypot(x - mean(x), y - mean(y)); % NOTE this centers on the actual distribution, allowing it to be not at detector centre
    if isnan(resolutions(ri).vertical_mean)
      resolutions(ri).vertical_mean = mean(x);
      resolutions(ri).vertical_std = std(x);
    else
      resolutions(ri).vertical_mean = mean([resolutions(ri).vertical_mean mean(x)]); % average the two masses (assumes this won't be run for more than two masses)
      resolutions(ri).vertical_std = mean([resolutions(ri).vertical_std std(x)]); % average the two masses (assumes this won't be run for more than two masses)
    end
    % TODO if optimize_potentials_objective.m adds a small penalty if the ions are not centered on the detector, include vertical_mean in addition to y_mean
    
    if ~isempty(x) && any(~isnan(x))
      if length(t) > 4 && length(r) > 4
        coeff_r(:,i) = ( (repmat(r,1,3) .^ repmat([0 2 4],length(r),1)) \ (t - t_mean_for_mass) )';
        coeff1_x(:,i) = ( (repmat(x,1,2) .^ repmat([0 1],length(x),1)) \ (t - t_mean_for_mass) )';
        ratio_wx_wy(1,i) = std(x) / std(y);
      else % too few hits to fit three coefficients
        coeff_r(:,i) = NaN;
        coeff1_x(:,i) = NaN;
        ratio_wx_wy(1,i) = NaN;
      end
      if show_figure && ~isempty(t)
        if ~show_x_instead_of_r
          % r instead of x, because azimuth-variation looks like degradiation of (x,t)- and (y,t)- but not (r,t)-plots
          x = r;
        end
        subplot(2,5,(row-1)*5 + 3);
        plot(x, t - t_mean_for_mass, style);
        hold on;
        h_lines(i) = plot(abscissa_y, polyval([coeff_r(end,i) 0 coeff_r(2,i) 0 coeff_r(1,i)], abscissa_y),  ['-' style(end)]); % quadratic & cubic fit
        plot(abscissa_y, polyval(coeff1_x([2 1],i), abscissa_y),  ['--' style(end)]); % linear fit
        ylabel('t - <t> / ns')
        if ~isempty(t_axis_range) && all(isfinite(t_axis_range)) && diff(t_axis_range)>0
          ylim(t_axis_range - t_mean_for_mass);
        end
        if resolutions(ri,1).kinetic_energy_used <= .2 %[eV]
          standard_det_halfrange = 5; %[mm]
        elseif resolutions(ri,1).kinetic_energy_used <= 1 %[eV]
          standard_det_halfrange = 10; %[mm]
        else
          %standard_det_halfrange = 30; %[mm]
          standard_det_halfrange = 22; %[mm]  detector is only 40 mm in diameter, thus show 44 rather than 60 mm diameter
        end
        
        if show_x_instead_of_r % (probably less useful for VMI, but to easier notice skewness in x distribution)
          if length(z_detector) == 1 % for old geometries without angled detector
            xlabel('vertical |x| / mm')
            if diff(minmax(x(:)')) < 2 * standard_det_halfrange
              xlim( mean(minmax(x(:)')) + [-1 1]*standard_det_halfrange );
            else
              if diff(minmax(x(:)')) < 6 * standard_det_halfrange
                xlim( mean(minmax(x(:)')) + [-3 3]*standard_det_halfrange );
              else % (old) show full range
                xlim([0 abscissa_y(end)])
              end
            end
          else % used for interesting geometries:
            xlabel('vertical x / mm')
            %if diff(minmax(x(:)')) < 2 * standard_det_halfrange
            %  xlim( mean(minmax(x(:)')) + [-1 1]*standard_det_halfrange )
            if max(abs(x(:))) <= 1 * standard_det_halfrange && ~isfinite(max_r)
              xlim( 0 + [-1 1]*standard_det_halfrange );
            else
              %if diff(minmax(x(:)')) < 6 * standard_det_halfrange
              %  xlim( mean(minmax(x(:)')) + [-3 3]*standard_det_halfrange )
              %else % (old) show full range
              %  xlim([abscissa_y(1) abscissa_y(end)])
              %end
              if isfinite(max_r)
                xlim([-max_r max_r]);
              else
                xlim([-22 22]); % as the high-energy range (44 mm diameter is a bit more than the detector)
              end
            end
          end
        else
          % When showing radius: only positive
          xlabel('refl.-radial sqrt(x^2 + y^2) / mm')
          xlim(xlim() .* [0 1]); 
          if max(x(:)) < standard_det_halfrange && ~isfinite(max_r)
            xlim([0 1]*standard_det_halfrange)
          else
            %if diff(minmax(x(:)')) < 3 * standard_det_halfrange
            %  xlim([0 3]*standard_det_halfrange)
            %else % (old) show full range
            %  xlim([0 abscissa_y(end)])
            %end
            if isfinite(max_r)
              xlim([0 max_r]);
            else
              xlim([0 22]); % detector diameter is 40 mm, show up to 44 mm diameter (22 mm radius)
            end
          end
        end
        
        i_of_same_mass = find(masses == flys(i).source.mass);
        if i == i2(end) || (~isempty(i1) && i == i1(end))
          % Last group of this mass is shown
          title(sprintf('Q_{yt}=%s, Q_{xt}=%s, \\sigmax/\\sigmay=%s', ...
            texformat_g10(mean(coeff1_y(2,i_of_same_mass)),2), ...
            texformat_g10(mean(coeff1_x(2,i_of_same_mass)),2), ...
            texformat_g10(mean(ratio_wx_wy(1,i_of_same_mass)),2) ));
        end
      end
      
    end % END x
    
    
    y0 = traj.y0(which) - nominal_y0;
    coeff_y0(:,i) = ( (repmat(y0,1,3) .^ repmat([0 2 4],length(y0),1)) \ (t - t_mean_for_mass) );
    % By analogy with Q_z for z_0, define a dimensionless Q_r0t such that dt/d|r0| = Q_r0t / v1
    % Initial case, assuming x=0 gave: v1 * dt = Q_r0t * dy_0 <==> scaled_dt = Q_r0t * dy_0
    %Q_r0t(i) = (abs(y0)*1E-3) \ scaled_dt; % least-squares fit of a deviation coefficient ("deblurring factor" for |y_0| to TOF)
    Q_r0t(i) = (hypot(y0, traj.x0(which))*1E-3) \ scaled_dt; % least-squares fit of a deviation coefficient ("deblurring factor" for sqrt(y_0^2+x_0^2) to TOF
    
    
    Q_y0t(i) = (y0*1E-3) \ scaled_dt; % least-squares fit of a deviation coefficient ("deblurring factor" for y_0 to TOF, relevant mainly in reflectron where central symmetry not guaranteed)
    if show_figure && ~isempty(t)
      subplot(2,5,(row-1)*5 + 4);
      %plot(traj.y0, traj.t - t_mean_for_mass, style);
      if length(z_detector) == 1
        plot(abs(y0), t - t_mean_for_mass, style); % dt[ns] = scaled_dt/v1 / 1E-9
      else
        plot(y0, t - t_mean_for_mass, style); % dt[ns] = scaled_dt/v1 / 1E-9
      end
      hold on;
      plot(abscissa_y0, polyval([coeff_y0(end,i) 0 coeff_y0(2,i) 0 coeff_y0(1,i)], abscissa_y0),  ['-' style(end)]);
      plot(abscissa_y0, Q_r0t(i) * abs(abscissa_y0)*1E-3 /v1/1E-9, ['--' style(end)]); % rescale curve according to dt[ns] = scaled_dt/v1 / 1E-9
      plot(abscissa_y0, Q_y0t(i) * abscissa_y0*1E-3 /v1/1E-9, [':' style(end)]); % rescale curve according to dt[ns] = scaled_dt/v1 / 1E-9
      ylabel('t - <t> / ns')
      xlabel('|y_0| / mm')
      if ~isempty(t_axis_range) && all(isfinite(t_axis_range)) && diff(t_axis_range) > 0
        ylim(t_axis_range - t_mean_for_mass);
      end
      if length(z_detector) == 1
        xlim([0 abscissa_y0(end)])
      else
        xlim([abscissa_y0(1) abscissa_y0(end)])
      end
      i_of_same_mass = find(masses == flys(i).source.mass);
      if i == i2(end) || (~isempty(i1) && i == i1(end))
        % Last group of this mass is shown
        %title(sprintf('<Q_{rt}>=%s', texformat_g10(mean(Q_r0t(i_of_same_mass)),4)));
        title(sprintf('Q_{rt}=%s, Q_{y0t}=%s', texformat_g10(mean(Q_r0t(i_of_same_mass)),4), texformat_g10(mean(Q_y0t(i_of_same_mass)),3)));
      end
    end


    dz0 = (traj.z0(which)-nominal_z0) * z_sign;
    % Using my Q_z "Wiley-McLaren deviation" definition: dt/dz_0 = Q_z / sqrt(2 q U_ext / m) = Q_z / v1
    % Q_z is dimensionless if all ingoing quantities are expressed in SI (base) units
    % Linear slope is Q_z. sqrt(2 q U_ext / m) * dt = Q_z * dz_0 <==> v1 * dt = Q_z * dz_0
    Q_z(i) = (dz0*1E-3) \ scaled_dt; % least-squares fit of a Wiley-McLaren deviation coefficient ("deblurring factor" for z_0 to TOF)
    if show_figure && ~isempty(t)
      subplot(2,5,(row-1)*5 + 5);
      plot(dz0, scaled_dt/1E-3, style);
      hold on;
      plot(full_dz0_range , Q_z(i) * full_dz0_range*1E-3 / 1E-3, ['--' style(end)]);
      xlabel('z_0 - <z_0> / mm')
      ylabel('Equivalent drift length offset: (t - <t>) \times v_1 / mm')
      xlim(full_dz0_range);
      if diff(t_axis_range) > 0
        ylim((t_axis_range - t_mean_for_mass)*1E-9 * v1 / 1E-3);
      end
      if i == i2(end) || (~isempty(i1) && i == i1(end))
        % Last group of this mass is shown
        title(sprintf('+Q_z=%s', texformat_g10(mean(Q_z(i_of_same_mass)),4)));
      end
    end
  end
  
  % TODO IMPROVEMENT: could add a plot with y(v_y0) i.e. to show velocity-map feasibility (by comparing to v_y0*t a magnification quantity can also be fitted & shown) 

  
  %if show_figure && all(isfinite(h_lines([i1 i2])))
  if show_figure && any(isfinite(h_lines([i1 i2])))
    labels = {};
    for i = [i2 i1]
      labels{i} = sprintf('%seV', texformat_SI(flys(i).source.energy));
    end

    if strcmp(labels(i1(1)),labels(i1(end))) && strcmp(labels(i2(1)),labels(i2(end)))
      % Seems same energy everywhere in diagram
      single_energy = [' ' labels{i1(1)}];

      % Monoatomic gases have 3 degrees of freedom, polyatomic also have rotational and vibrational,
      % so if one knew their total energy one should divide by more to get the temperature.
      % However what flys(i).source.energy contains just the TRANSLATIONAL kinetic energy,
      % so maybe division by 3 is still correct?
      % https://en.wikipedia.org/wiki/Degrees_of_freedom_(physics_and_chemistry)
      % https://en.wikipedia.org/wiki/Kinetic_theory_of_gases
      % https://en.wikipedia.org/wiki/Heat_capacity_ratio
      translational_degrees_of_freedom = 3; 
      temperature = 2 * flys(i).source.energy*eV / (translational_degrees_of_freedom * k_Boltzmann);
      if temperature >= 999.5 %[K]
        single_energy = sprintf('%s (%.0f kK)', single_energy, temperature/1E3);
      else
        single_energy = sprintf('%s (%.0f K)', single_energy, temperature);
      end
    else
      single_energy = '';
      
      subplot(2,5,1);
      h = legend(h_lines(i2), labels(i2), 'Location','SE');
      pos = get(h,'Position');
      pos(1) = 0.2; % move legend over to the TOF-histogram diagram (leftmost)
      set(h,'Position', pos, 'Color','none'); 

      subplot(2,5,5);
      h = legend(h_lines(i1), labels(i1), 'Location','SE');
      pos = get(h,'Position');
      pos(1) = 0.2; % move legend over to the TOF-histogram diagram (leftmost)
      set(h,'Position', pos, 'Color','none'); 
    end

    handles.central_annotation = annotation('textbox', [0.28 0.48 0.6 0.04], 'String', ...
        strrep(sprintf('By spacing%s: %.1fu with %seV %s', single_energy, resolutions(ri).mass_extrapolated_by_spacing, ...
                texformat_SI(mean(resolutions(ri).kinetic_energy_used)), annotation_suffix), ...
               '\mu','µ'), ...
        'EdgeColor','none', 'Interpreter','none');

  %   coeff_y
  %   coeff_y0
  %   Q_r0t
  %   Q_z
  %   coeff1_x

    if isfield(flys(1), 'workbench_name')
      n = [' ' flys(1).workbench_name];
    else
      n = '';
    end
    if size(flys(1).rays,1) == 1 && flys(1).imaging_dimensions == 3
      % (hints that) Isotropic velocity distribution, not constrained to a plane (used for mass resolution when no eVMI)
      % Let file name start with an "i" instead of "m".
      type = 'i';
    else
      type = 'm';
    end
    if show_x_instead_of_r
      type = [type 'x'];
    end

    if flys(1).potentials(1) == 0 && length(flys(1).potentials(2:end)) < 18
      % if first potential is zero and length is below 18 % (i.e. not "both" geometry): skip first potential
      figure_filename = [type n ' -' strrep(strrep(sprintf(' %+.0f', round(flys(1).potentials(2:end)')), '-',''), '+0','0')]; % skip (1) which is the grounded walls
    else % show all
      figure_filename = [type n ' -' strrep(strrep(sprintf(' %+.0f', round(flys(1).potentials')), '-',''), '+0','0')]; % skip (1) which is the grounded walls
    end
    
    if all(resolutions(ri).mass_used ~= 200) || any(resolutions(ri).kinetic_energy_used ~= 1)
      %figure_filename = sprintf('%s %.3geV@%.0fu', figure_filename, resolutions(ri).kinetic_energy_used(1), resolutions(ri).mass_used(1)); % energy at end of name
      figure_filename = sprintf('%s %seV@%.0fu%s', type, strrep(texformat_SI(resolutions(ri).kinetic_energy_used(1),2,false), ' ',''), resolutions(ri).mass_used(1), figure_filename(2:end)); % energy at beginning of name
      % Shorten by using mV and ku 
      figure_filename = strrep(figure_filename, '000u', 'ku'); % [kDa] is more SI, but the aim is to shorten the file name ;-)
    end
    if ~isempty(source_FWHM_y_str)
      figure_filename = [figure_filename ' ' source_FWHM_y_str];
    end
    figure_filename = strrep(figure_filename, 'both8eu_', ''); % shorten long filenames when 19 electrodes in both8eu_B2

    % Adjust axis positions, default doesn't use the full figure width
    h_axes = get(8, 'Children');
    p = cell2mat(get(h_axes, 'Position'));
    for i = 1:length(h_axes)
      row = p(i,2) > 0.5;
      col = floor(p(i,1) / 0.16);
      set(h_axes(i), 'Position', [0.055+col*0.195, 0.058+0.51*row, 0.15, 0.37]);
    end
    
    export_fig_multiformat(8, figure_filename, 'PNG');
  end
  
  % End of analysis for one [i2 and possibly i1] pair.
  % Mark these masses as processed so next iteration will use pair of lower mases, if any.
  masses(masses >= mass2) = -mass2;
  if ~isempty(mass1)
    masses(masses >= mass1) = -mass1;
  end
end
  
geometric_coefficients = struct('coeff_y',coeff_y, 'coeff_y0',coeff_y0, 'coeff_r',coeff_r, ...
  'coeff1_y',coeff1_y, 'coeff1_x',coeff1_x, 'Q_r0t', Q_r0t, 'Q_z',Q_z, 'ratio_wx_wy',ratio_wx_wy);


