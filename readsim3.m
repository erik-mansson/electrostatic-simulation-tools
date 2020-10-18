% flys = readsim3()
% fly  = readsim3(position)
% fly  = readsim3(position, filename)
%
% Reads one or several entries (groups) from a record file written by SIMION.
% Configure SIMION's Data recording by loading record3.rec
% and in the bottom right box setting the output file to
% the (full) path of the filename that will be read here.
%
% For this version, readsim3, the Lua workbench script is able to run different 
% kinetic energies in one "Fly'm" command. This means that raw data from all entries 
% is given in the beginning (or not given at all) and then followed by the summaries (rays)
% of each entry. Voltage info is only given once, at the beginning.
% (In the previous version, readsim2 (used together with runsim2 & commandline.lua),
% more text was repeated for each entry (e.g. voltages) and raw data (individual 
% particle coordinates) were given within each entry.)
%
% NOTE: the next SIMION run appends to the file, so it eventually
% needs to be deleted to clear old output! If SIMION is launched
% via runsim3 any old output is first deleted.
% 
%
% PARAMETERS
%   position   (Default: 0)
%              A text file may contain multiple SIMION-runs.
%              position=0 means that we read read all of the runs, returned as a struct array.
%              A nonzero index (e.g. 1), or an array of indices, means that those specific runs are read.
%   filename   The file to read. Default 'sim.txt' but normally you'd want to give an absolute path to it.
%
% RETURN
%   A "fly"-structure or a column vector of such, containing the following:
%   fly.U_MCP: Potential on detector front [V]
%   fly.U_last: Potential on electrode closest to detector [V]
%   fly.U_1   Potential on electrode most distant from detector [V]
%             U_1 is redundant, because U_1 = U_last + sum(fly.V)
%   fly.V:    (electrode_count-1)-by-1 voltage array between electrodes. [V]
%   fly.source.ion_count
%   fly.source.energy      [eV]
%   fly.source.mass        [u]=[a.m.u.]
%   fly.source.charge      [q]=[elementary charge]
%   fly.source.point_count Number of different source points in space
%   fly.source.pattern:    How are the source points spread?
%                            0: along a radial line (at constant lab-z, i.e. constant SIMION-x)
%                            1: a pattern with both radial and axial (lab-z i.e. SIMION-x) coordinate varying.
%   fly.source.spacing:    Related to the spacing between adjacent source
%                          points, the exact meaning depends on fly.source.pattern.
%   fly.rays:              Matrix (direction_count-by-8or10) with info about the hits
%                          For instance, strcmp(fly.ray_columns,'y_mean'), tells which column that contains the y-mean.
%   fly.rays(direction_index,:) = [theta, hit_count, y_mean, y_std, y_range, ...
%                                  t_mean, t_std, hit_ratio, y2_std, y2_std_rel]
%                          or   = [theta, hit_count, y_mean, y_std, y_range, ...
%                                  t_mean, t_std, hit_ratio]
%     The distance (y) unit is [mm], time (t) unit is [ns] and the theta unit is degrees [pi/180].
%     The relative quantity y2_std_rel is a fraction in the range 0 to 1, not a percentage. Its inverse is the y^2-resolving power.
%   fly.misses             The number of ions that missed the detector
%
%   fly.overall            A structure with overall resolution info, calculated in a
%                          rather rough way in Simion (lens_manual2.sl).
%   fly.overall.y_std      2-norm average of the fly.rays.y_std (equal weight for each ray)
%                          where a minimum of 0.5mm is used for each fly.rays.y_std.
%                 formula: norm(max(rays(:,4),0.5))/sqrt(size(rays,1))
%   fly.overall.y_rel      (approximative) average of 2-norm average for the 
%                          relative error in the y-direction (y_std/y_mean).
%                          A minimum value for y_mean is used to accept y=0
%    approximate formula?: norm(max(rays(:,4),0.5)./max(abs(rays(:,3)),7))/sqrt(size(rays,1))
%   fly.overall.t_std      2-norm average of the fly.rays.t_std (equal weight for each ray)
%                          with some minimum limits...
%   fly.overall.t_rel      (approximative) 2-norm average for the 
%                          relative error in the time-direction (t_std/(t_mean-t_0))
%                          with some approximations and minimum limits...
%                          Some directions are excluded from this calculation (e.g. theta=90)
%   fly.overall.tot_rel    the 2-norm average of y_rel and t_rel
%                 formula: norm([y_rel t_rel])/sqrt(2)
%
%   fly.trajectories       Struct with N-by-1 fields describing individual particles, or empty if no such info was contained in the file.
%                          For runsim3, the varaible-value pair 'trajectories',true can be given to ensure trajectory info is logged.
%                          NOTE: the x, y, z convention is swapped with respect to the SIMION coordinates (Z,Y,X).
%                          SIMION has X as the time-of-flight axis while we prefer to call it z here, 
%                          consequently the vertical axis is here called x. In both coordinate systems, y (and Y) is the optical axis.
%                            y [mm]  final transverse position (along optical axis)
%                            t [ns]  time of flight
%                            mass [u]
%                            charge [e]
%                            energies [eV]       initial kinetic energy
%                            elevation [degrees] of initial velocity
%                            z0 [mm] initial axial coordinate (absolute) (time of flight axis)
%                            y0 [mm] initial transverse position (absolute)
%                            z  [mm] final axial coordinate (absolute) (time of flight axis) -- for telling whether
%                                    particle hit detector or something else (although this axis is X in SIMION).
%                            x [mm]  final vertical transverse position (although this axis is Z in SIMION)
%                            x0 [mm] initial vertical transverse position (although this axis is Z in SIMION)
%   (only for the last fly):
%   flys(end).geometry      .detector_A_coordinates [mm]
%                                                     The [z-centre, y-centre, x-centre] for the VMI detector.
%                                                     Note that the TOF-axis is first and vertical last (as the SIMION-coordinates are ordered).
%                           .detectorB_xSIMION_zMatlab [mm] coordinate for centre of ion (reflectron) detector along 
%                                                     time of flight axis, also called detectorB_x or MCP_centre_x
%                           .detectorB_linecoeff [mm] also called detectorB_x in Lua and reflectron optimization programs.
%                                                     Defined as [k, -1, -d, MCP_centre_y]
%                                                     where k=-1/tan(det.angle=6deg) 
%                                                     and d = k*MCP_centre_x - 1*MCP_centre_y
%
% EXAMPLE: reads the last run in the default filename
%   fly=readsim2(-1); 
%   plot(fly.trajectories.t, fly.trajectories.y, '.')
%   fly(end).V 
% EXAMPLE: show the coloured TOF,|y| diagram for the last run
%   fly=readsim2(-1); show(fly);
% EXAMPLE: show the TOF,y diagram for all runs in the file
%   flys=readsim2; show(flys);
% EXAMPLE: print one column of essential info (assuming V1..V6 are identical)
%   % the U_1 is in [kV] and the last value(s) are the "resistor ratio(s)"
%   fly=readsim(-1); m=min(abs(fly.V(abs(fly.V)>0)));
%   [fly.U_last;fly.U_1/1e3;mean(fly.V(1:6));fly.V(7);fly.V(8:end);fly.overall.y_std;fly.overall.t_std;
%    fly.overall.tot_rel*100;fly.overall.y_rel*100;fly.overall.t_rel*100;fly.misses;unique(abs(fly.V(abs(fly.V)>m)))/m]
% EXAMPPLE: mean total relative error for all runs in a file
%   flys=readsim; mean(collect_field(collect_field(flys,'overall'),'tot_rel'))
%
% SEE ALSO
%  runsim3, script_VMI_sim_and_calib.m, optimize_potentials_objective.m
function flys = readsim2(position_within_file, filename)
if nargin < 1
  position_within_file = 0; % get all
end
if nargin < 2
  filename = 'sim.txt'; % default for where to find output log file
end
flys = NaN; % in case of error

if ~exist(filename,'file')
  error('File missing: %s', filename); return
end
constants;
% Read entire file into a string
fp = fopen(filename, 'r');
file_contents = char(fread(fp,'char')');
fclose(fp);

% Skip past any comments (info not to be parsed by this program, those lines start by "# ")

% Find start positions of the flyings in the file, by searching for texts
run_starts = strfind(file_contents, 'U_AMCP'); % Occurs near beginning of each SIMION "Fly'm".
summary_starts = strfind(file_contents, sprintf('\nFlew')) + 1; % From version 3 of the workbench script, there may be multiple groups (each writing a "Flew"-summary) within one SIMION "Fly'm")
run_index = NaN(size(summary_starts)); % for each summary, tell which run it belongs to
for r = 1:length(run_starts)
  run_index(summary_starts >= run_starts(r)) = r;
end
if any(isnan(run_index)) || any(diff(run_index) < 0) || isempty(summary_starts)
  error('Failed to parse high-level structure in %s. Maybe no Lua-script output? This can happen because segment.terminate() is not run if the last particle leaves the simulation volume without hitting anything.', filename);
end
% NOTE: could let output have a struct of runs too (holding voltages and raw data without grouping by group/summary, with flys referring to them by index, but it would break some compatibility with readsim2.
% Now, emulating runsim2-like output by copying voltage data to each group/summary.
% IMPROVEMENT: read raw data and try to split it by group to be runsim2-compatible.

% Negative positions means to count from the end (-1 ==> last run)
position_within_file(position_within_file < 0) = length(summary_starts) + 1 + position_within_file(position_within_file < 0);
% A zero anywhere means to use all runs
use_all = any(position_within_file == 0);
returned_runs = 0;
previous_run = 0;

for index = 1:length(summary_starts)
  % Select the part(s) to use:
  if ~use_all && ~any(position_within_file == index)
    continue; % Do not use this summary
    % NOTE: not processing all runs & summaries may give problem with the (new for readsim3) ion_recording treatment.
  end
  counts_per_ray = []; %(default, if same number of particles in each ray/bunch/direction)

  if run_index(index) ~= previous_run
    % Need to parse and read electric potentials for this run. (NOTE: they are kept in simple variables,
    % not in a struct for indexed recall later. This is sufficient because run_index is monotonously non-decreasing.)
    previous_run = run_index(index);
    
    % Default in case not found:
    electrodes = struct('A_MCP',[], 'A_LAST',[], 'A_FIRST',[], 'A_FREEGROUP_COUNT', [], 'B_MCP',[], 'B_LAST',[], 'B_FIRST',[], 'B_FREEGROUP_COUNT',[], 'WALLS',[]);
    trajectories_ungrouped = struct('y',[], 't', [], 'mass',[], 'charge',[], ...
                                    'energies',[], 'elevation',[], 'z0', [], 'y0', [], 'z', [], 'x', [] );

    % The potentials are given before the first summary of the run:
    potential_string = file_contents(run_starts(run_index(index)):(summary_starts(find(run_index == run_index(index), 1))-1));
    
    [values,count,~,nextindex] = sscanf(potential_string, 'U_AMCP = %g, U_Alast = %g, U_BMCP = %g, U_Blast = %g, A_free=%g, A_plateau=%g\n', 6);
    % Note that when using DoNotSetPotentials=-1, the A_free, A_plateau and B_FREEGROUP_COUNT will be NaN.
    if count == 6
      U_AMCP = values(1); U_Alast = values(2);
      U_BMCP = values(3); U_Blast = values(4);
      A_free = values(5); A_plateau = values(6); 
      potential_string = potential_string(nextindex:end);
      [V,count,~,nextindex] = sscanf(potential_string, ' U[%d] = %g, V = %g\n');
      if count >= 3 && mod(count, 3) == 0 % if an integer number of rows were parsed
        % OK
        potential_index = V(1:3:end); % gives the "historical" indexing of U[] in sim.txt (the line "U[8] = -190,  V = -10" means potentials(potential_index==8) = -190).
        potentials = V(2:3:end); % potentials, where 
        V = V(3:3:end); % voltages (differences between adjacent potentials). NOTE: V(end) is NaN, because undefined to which electrode a difference should be taken (could use chamber walls though...)
        if potential_index(end) ~= length(potential_index) || any(potential_index ~= [1:length(potential_index)]')
          % From version 3 requiring to get full list, starting at index 1.
          error('Potentials of some electrodes are missing. Got indices %s from\%s', mat2str(potential_index'), potential_string);
        end
        potential_string = potential_string(nextindex:end);
      else
        error('Unrecognized format for voltages in run %d:\n%s', index, potential_string); return
      end
      clear potential_index % since all indices from 1 and up are present, this array is not needed
      
      [I,count,~,nextindex] = sscanf(potential_string, 'Indices: A_MCP=%d, A_LAST=%d, A_FIRST=%d, B_MCP=%d, B_LAST=%d, B_FIRST=%d, B_FREEGROUP_COUNT=%g, WALLS=%d\n', 8);
      if count == 8
        % Explains the meaning of some indices in the potential array, e.g. potential(electrodes.A_FIRST)
        electrodes = struct('A_MCP',I(1), 'A_LAST',I(2), 'A_FIRST',I(3), 'A_FREEGROUP_COUNT', A_free, 'B_MCP',I(4), 'B_LAST',I(5), 'B_FIRST',I(6), 'B_FREEGROUP_COUNT',I(7), 'WALLS',I(8));
        potential_string = potential_string(nextindex:end);
      else
        warning('Unrecognized format for electrode indices.');
        % (Using the default electrodes set above)
      end

      % NOTE IMPROVEMENT: put in Lua, this formula assumes the needle is halfway between B_first and A_first.
      U_ext = (potentials(electrodes.B_FIRST) - potentials(electrodes.A_FIRST)) / 2;
      W_ext = abs(U_ext);
      % NOTE: at this stage the particle charge is not known. U_ext, U_tube and U_tweak will be updated later (see below).
      V_Abend = (potentials(electrodes.A_FIRST) - potentials(electrodes.A_LAST));
      V_Atweak = (potentials(electrodes.A_FIRST+max([A_free 0])) - potentials(electrodes.A_FIRST));
      if electrodes.B_FIRST == 0 || electrodes.B_LAST == 0 % to allow small geometries where no first/end electrode on B-side
        V_Bbend = NaN;
      else
        V_Bbend = (potentials(electrodes.B_FIRST) - potentials(electrodes.B_LAST));
      end
      V_Btweak = (potentials(electrodes.B_FIRST-max([electrodes.B_FREEGROUP_COUNT 0])) - potentials(electrodes.B_FIRST));
      
    else
      % No info about potentials, read only basic particle info
      U_AMCP = NaN; U_Alast = NaN;
      U_BMCP = NaN; U_Blast = NaN;
      A_free = NaN; A_plateau = NaN; 
      potential_string = potential_string(nextindex:end);
      V = [];
      potentials = [];
      U_ext = NaN; W_ext = NaN;
      V_Abend = NaN;
      V_Atweak = NaN;
      V_Bbend = NaN;
      V_Btweak = NaN;
      warning('Unrecognized format for potentials in run %d (for summary %d).', run_index(index), index);
%       error('Unrecognized format for log file for summary %d (in run %d):\n%s', index, run_index, potential_string);
    end
    
    % Read lens parameters
    % The symbol U is used for "signed potential energy difference" ("work with sign"), which is positive when accelerating the particel along the spectrometer axis (in a particular direction). For U_tweak positive means (focusing, demagnifying) repelling wrt. extractor.
    % The symbol V is used for "absolute electric potential" ("voltage with respect to ground")
    % NOTE: at this stage the particle charge is not known. U_ext, U_tube and U_tweak will be updated later (see below, per summary).
    % Now using signed W_ext, so that extraction work will be positive only for the side (A or B) to which particle of given sign is actually extracted!
    lens_parameters_A = struct('free',A_free, 'plateau',A_plateau, 'U_ext',U_ext, 'U_bend',V_Abend, 'bend',V_Abend/U_ext, 'U_tweak',V_Atweak, 'tweak1',V_Atweak/U_ext);
    [lens_parameters_A, potential_string] = read_lens_parameters(potential_string, lens_parameters_A);
    lens_parameters_A.V_last = U_Alast; lens_parameters_A.V_MCP = U_AMCP; % absolute potentials (voltage wrt. ground, independent of particle charge)
    lens_parameters_B = struct('free',1,      'plateau',0,         'U_ext',-U_ext, 'U_bend',V_Bbend, 'bend',V_Bbend/U_ext, 'U_tweak',V_Btweak, 'tweak1',V_Btweak/U_ext);
    [lens_parameters_B, potential_string] = read_lens_parameters(potential_string, lens_parameters_B);
    lens_parameters_B.V_last = U_Blast; lens_parameters_B.V_MCP = U_BMCP; % absolute potentials (voltage wrt. ground, independent of particle charge)

    if isnan(lens_parameters_A.E_extraction)
      % There was no info about parameters (this should only happen if DoNotSetPotentials was nonzero, otherwise an error).
      if isnan(A_free)
        % Seems like DoNotSetPotentials == -1, i.e. the potential array was set directly rather than via the few parameters
        % Compute some things that read_lens_parameters will have failed, but that are possible to determine:
        lens_parameters_A.extractor_repeller_ratio = (potentials(electrodes.A_FIRST)-potentials(electrodes.A_LAST)) / (potentials(electrodes.B_FIRST)-potentials(electrodes.A_LAST)); % VMI style extractor/repeller ratio (with respect to detector&tube, ignoring that lens may deviate)
        if electrodes.B_FIRST == 0 || electrodes.B_LAST == 0 % to allow small geometries where no first/end electrode on B-side
          lens_parameters_B.extractor_repeller_ratio = NaN;
        else
          lens_parameters_B.extractor_repeller_ratio = (potentials(electrodes.B_FIRST)-potentials(electrodes.B_LAST)) / (potentials(electrodes.A_FIRST)-potentials(electrodes.B_LAST));
        end
      else
        warning('No lens parameter summary was read.')
      end
    end

    imaging_dimensions = 3; % the only mode in readsim2
    counts_per_ray = sscanf(potential_string, '# Using VMI-style 3D-isotropic, both p_Z signs within each ray. Count per ray: [%[0-9, ]]\n', 1);
    if ~isempty(counts_per_ray)
      % Using _isotropic = 3 or 4, where the number of particles flewn per ray/bunch is not constant.
      counts_per_ray = transpose(str2num(counts_per_ray));
      % 45 and 135 degree elevation are handled as the same ray/bunch
      imaging_dimensions = 2;      
    elseif ~isempty(strfind(potential_string, '# Using VMI-style distribution')) % possible (but not required) in readsim3
      % Then 45 and 135 degree elevation are handled as the same ray/bunch
      imaging_dimensions = 2;
    end

    % Skip past one (or zero) block of comments or format descriptions
    % (a sequence of lines that are empty or start with '"' or '# ').
    potential_string = strtrim(regexprep(potential_string, '\n*((# |")[^\n]+\n|\r?\n)+', '', 1));
    
    if isempty(potential_string)
      ion_recording = [];
    else
      % IMPROVEMENT: read raw data (coordinates) too
      % Read ion data recording into matrix (when recording_enable=1)
      [ion_recording,count,~,nextindex] = sscanf(potential_string, '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n', Inf);
      if mod(count, 11) == 0 && count >= 11
        % Recording OK, by using record3.rec, including third coordinate (z in SIMION)
        potential_string = potential_string(nextindex:end);
        ion_recording = reshape(ion_recording,11,length(ion_recording)/11)';
        % TODO for the readsim3 type of log, where multiple summaries (groups) within the same fly: use the summaries later to split the ion_recording into part per group. E.g. simply use the number of particles and summaries!
      else
        warning('readsim3:trajectories', 'The recording is not by record3.rec format, or there is some special value like "-1.#IND" indicating computation error, e.g. when defining velocity. Will try reading older formats, just in case.');
        
        [ion_recording,count,~,nextindex] = sscanf(potential_string, '%g,%g,%g,%g,%g,%g,%g,%g,%g\n', Inf);
        if mod(count, 9) == 0 && count >= 9
          % Recording OK, by using record.rec, without SIMION z-coordinate and velocity
          potential_string = potential_string(nextindex:end);
          ion_recording = reshape(ion_recording,9,length(ion_recording)/9)';
          % TODO for the readsim3 type of log, where multiple summaries (groups) within the same fly: use the summaryies later to split the ion_recording into part per group. E.g. simply use the number of particles and summaries!
        else % Failed to read recorded ion data
          count
          potential_string(1:min([200 length(potential_string)]))
          ion_recording
          count = 0
          ion_recording = [];
          warning('Failed to read raw coordinate data.');
        end
      end
    end
  end
  clear potential_string
  
  % Read a summary
  if index < length(summary_starts)
    s = file_contents(summary_starts(index):(summary_starts(index+1)-1)); % the relevant string part of the file
  else
    s = file_contents(summary_starts(index):end);
  end
  
  % Read basic info about the ion source
  [values,count,~,nextindex] = sscanf(s, 'Flew %d ions @ %g eV, %g u, %g q. %g points: pattern %d, spacing %g mm\n', 7);
  if count ~= 7
    [values,count,~,nextindex] = sscanf(s, 'Flew %d ions @ %g eV, %g u, %g q. %g points: gaussian std [%f,%f,%f] mm\n', 8);
    if count ~= 8
      error('Unrecognized format for hit info per direction for summary %d (in run %d):\n%s', index, run_index(index), s(1:min([200 length(s)])));
    end
  end
  charge = values(4);
  source = struct('ion_count',values(1), 'energy', values(2), ...
    'mass', values(3), 'charge', charge, ...
    'point_count',values(5), 'pattern',values(6), 'spacing',values(7));
  if count == 8
    % The Gaussian-in-FLY2 variant
    source.pattern = -1;
    source.spacing = values(6:8)'; % use a 1-by-3 vector of standard deviations instead of a spacing
  end
  s = s(nextindex:end); % Go to position after source info
  direction_count = ceil(source.ion_count / source.point_count);
  
  % Read info about the ray (trajectory bunch) for each direction.
  % Columns: rays(direction_index,:) = [theta, hit_count, y_mean, y_std, y_range, t_mean, t_std, hit_ratio, y2_mean, y2_std_rel]. Distances in [mm], times in [ns].
  columns = {'theta', 'hit_count', 'y_mean', 'y_std', 'y_range', 't_mean', 't_std', 'hit_ratio', 'y2_std', 'y2_std_rel'};
  % NOTE: the last two columns (for y^2) were not present in readsim1 & readsim2, reading old files will thus fail.
  % %*[s 0-9.%] skips past the column which has relative TOF-separation for some, but not all, directions (in 3D mode).
  %    "s" is included to capture the "s" in "ns" and avoid failure for the direction without any value (matching empty strings is not supported).
  if imaging_dimensions == 2 %(_isotropic = 2, 3 or 4)
    [values,count,~,nextindex] = sscanf(s, '%g d&: %d ions; y %g std %g, range %g mm; t%g std%g n%*[s 0-9.%$u]; y2 std%g mm^2%g%%%*[ \r\n]', 9 * direction_count);
  else % default: 3D (_isotropic=0 or 1)
    [values,count,~,nextindex] = sscanf(s, '%g d.: %d ions; y %g std %g, range %g mm; t%g std%g n%*[s 0-9.%$u]; y2 std%g mm^2%g%%%*[ \r\n]', 9 * direction_count);
  end
  if mod(count, 9) == 0 && count > 0
    rays = reshape(values,9,direction_count)';
    rays(rays(:,2)==0, 3:end) = NaN; % Use NaN for rays without any detector hit
    % Move the new y^2-info to index 9 and 10 to let 1 to 8 be compatible with readsim2.
    rays(:,9:10) = rays(:,8:9);
    rays(:,10)  = rays(:,10) / 100; % Convert from percentage (0-100) to fraction (0-1). Computed using y2_limited = max([MIN_R_REL^2 det_y2_mean(ion_group)(dir)]) to avoid dividing by too small y^2-value.
    if isempty(counts_per_ray)
      rays(:,8) = rays(:,2)/source.point_count; % create hit_ratio=hit_count/point_count (should equal 1)
    else
      if length(counts_per_ray) ~= size(rays,1)
        error('Got particle counts for %d rays, when there are %d rays.', length(counts_per_ray), size(rays,1));
      elseif sum(counts_per_ray) ~= source.ion_count
        error('Got particle counts that sum to %d while %d particles have been flown.', sum(counts_per_ray), source.ion_count);
      end
      rays(:,8) = rays(:,2) ./ counts_per_ray; % create hit_ratio=hit_count/point_count (should equal 1)
    end
    s = s(nextindex:end); % skip past the values
    %nextindex = strfind(s, 'Tot sy'); nextindex = nextindex(1); s = s(nextindex:end); % skip past the space and last newline, to find summary
  else
    error('Unrecognized format for direction (ray) in summary %d (in run %d):\n%s\n', index, run_index(index), s(1:min([200 length(s)])));
  end
  
  % Read overall resolution info:
  [values,count,~,nextindex] = sscanf(s, '%*[ \r\n]Tot sy %g mm, <sy/y> %g%%. st %g ns, <st/t> %g%%. <1/R> %g%%', 5);
  overall = struct('y_std', values(1), ...
    'y_rel', values(2)/100, ... % Convert from percentage (0-100) to fraction (0-1). Computed using y_limited = max([MIN_R_REL abs(det_y_mean(ion_group)(dir)]) to avoid dividing by too small |y|-value.
    't_std', values(3), 't_rel', values(4)/100, ...
    'tot_rel', values(5)/100); % Convert from percentage (0-100) to fraction (0-1).
  s = s((nextindex+1):end); % skip past the values. +1 since the trailing % sign seems to not be accounted for by the previous matching.
  [values,count,~,nextindex] = sscanf(s, '%*[ \r\n]Tot s(y^2)%g mm^2, <s(y^2)/y^2>%g%%', 2);
  if count == 2
    % Now y^2 statistics is present too. Since "radial kinetic energy" is proportional to y^2, this should give a more useful estimate of relative kinetic energy resolution.
    overall.y2_std_mean = values(1); % This is mean(std(y^2)), not root(mean(variance(y^2))). Thus a difference from .y_std and .t_std
    overall.y2_rel = values(2) / 100; % Convert from percentage (0-100) to fraction (0-1). Computed using y2_limited = max([MIN_R_REL^2 det_y2_mean(ion_group)(dir)]) to avoid dividing by too small y^2-value.
  else % Unexpectedly missing the y^2 info
    overall.y2_std_mean = NaN;
    overall.y2_rel = NaN;
  end
  
  if length(s) >= nextindex && s(nextindex)=='%'
    s = s(nextindex:end); % skip past the values and trailing % sign
  else
    s = s(nextindex:end); % skip past the values
  end
  % If s becomes empty it just means there was nothing extra logged before next summary_starts (that and following text has not been put in the variable s to begin with)
  di = strfind(s, 'Detector coordinates:');
  if ~isempty(di)
     s = s(di:end); % skip past any other stuff in the residual string (e.g. "Approximative max mass resolved")
  end
  [values,count,~,nextindex] = sscanf(s, 'Detector coordinates: A %f, %f, %f; B %f, %f, %f, %f, %f;', 8);
  if count == 8
    % New lua scripts (2020-02-25) print some geometry-specific info at the end (after last group)
    geometry = struct();
    geometry.detector_A_coordinates = transpose(values(1:3)); % x,y,z in SIMION: z,y,x in Matlab; (TOF axis, optical axis, vertical): detectorA_x, detectorA_y, VERTICAL_CENTRE in Lua
    geometry.detectorB_xSIMION_zMatlab = values(4); % x in SIMOIN, z in Matlab, "TOF axis": detectorB_x in Lua
    geometry.detectorB_linecoeff = transpose(values(5:end)); % Lua detectorB_linecoeff described as = [k, -1, -d, MCP_centre_y]: k=-1/tan(det.angle=6deg), d = k*MCP_centre_x - 1*MCP_centre_y;
  else
    geometry = []; % use an empty array so that isempty(geometry) will return true (it returns false for an empty struct)
  end
  
  returned_runs = returned_runs + 1;

  if ~isempty(ion_recording)
    % Copied from show.m: extract start and hit info per individual
    % trajectory (two adjacent rows from ion_recording):
    if size(ion_recording,2) == 9
      % If the old recording.rec was used, expand by inserting NaN as x (SIMION z) and v_z (SIMION v_z) coordinates
      ion_recording = [ion_recording(:,1:5) NaN(size(ion_recording,1),1) ion_recording(:,6:9) NaN(size(ion_recording,1),1)];
      has_x = false;
    else
      has_x = true;
    end
    
    start_indexes = find(ion_recording(:,1) == 0);
    starts = ion_recording(start_indexes, :);
    hits = ion_recording([start_indexes(2:end)-1;end], :);
    TOF = hits(:,1) * 1000; %convert to [ns]
    ion_mass = starts(:,2);
    
    % The number of digits printed gives a quite rough value for electron mass. Replace it with the more accurate value (may deviate slightly due to different CODATA-versions in contants.m and SIMION).
    electron_mass_in_u = electron_mass/atomic_mass_unit;
    electron_like = starts(:,3)==-1 & abs(ion_mass-electron_mass_in_u)<1E-10;
    ion_mass(electron_like) = electron_mass_in_u;
    
    vt = starts(:,7); %[mm/us]=[1e3 m/s] initial speed
    %Simion doesn't use all 360 degrees for elevation: elev = starts(:,8);
    if has_x
      elev = rad2deg(atan2(hypot(starts(:,11),starts(:,10)), starts(:,9)));
      azimuth = rad2deg(atan2(starts(:,11), starts(:,10)));
    else
      elev = rad2deg(atan2(starts(:,10),starts(:,9)));
      azimuth = NaN(size(starts,1),1);
    end
    % One entry per ion trajectory, gathering hit and start data
    % NOTE: X and Z are swapped from SIMION coordinates to the readout names.
    % In SIMION: x is the TOF axis but it is called trajectories.z here,
    % In SIMION: z is the vertical axis (polarization) but it is called trajectories.x here.
    % In SIMION: y is the optical axis, it is called trajectories.y.
    % NOTE: when plotting VMI image with optical axis horizontal, it means y is horizontal and x vertical...
    % NOTE: In the SIMION GUI's third choice of projection, the vertical axis (SIMION z, trajectories.x) is shown increasing downwards (but also the legend axis points that way).
    %       This is just a surprising visualization choice. The coordinate system is still as expected, and my both8... models (and other) are built so that positive vertical is up in reality.
    trajectories_ungrouped = struct('y',hits(:,5), 't', TOF, 'mass',ion_mass, 'charge',starts(:,3), ...
      'energies',atomic_mass_unit * ion_mass .* (vt*1e3).^2 / 2 / eV, ...
      'elevation',elev, 'z0', starts(:,4), 'y0', starts(:,5), 'z', hits(:,4), ...
      'x', hits(:,6), 'x0', starts(:,6), 'azimuth',azimuth );
    
    % Reset the ion_recording, so the same data won't be repeated in every group of this run.
    % NOTE: For runsim2 there could be only one group per run but for runsim3, the ion_recording
    % is printed once although it contains particles from all the groups of the run.
    % Thus a de-mixing procedure will be performed after this if-staement, moving some from trajectories_ungrouped to each group.
    ion_recording = [];
  end
  if isempty(trajectories_ungrouped.t)
    trajectories = struct('y',[], 't', [], 'mass',[], 'charge',[], ...
      'energies',[], 'elevation',[], 'z0', [], 'y0', [], 'z', [], 'x',[], 'x0',[], 'azimuth',[] );
  else
    
    % NOTE: For runsim2 there could be only one group per run but for runsim3, the ion_recording
    % is printed once although it contains particles from all the groups of the run.
    % Thus a de-mixing procedure is performed here, extracting particle trajectory coordinates belonging to current group.
    
    % This direct search works only when initial_vy=initial_v_vertical=0, i.e. the initial energies are as expected:
%     which_belong_here = trajectories_ungrouped.charge == source.charge ...
%                       & abs(trajectories_ungrouped.mass - source.mass) < 5E-5 ...
%                       & abs((trajectories_ungrouped.energies+0.05) ./ (source.energy+0.05) - 1) < 1E-4;
    % Alternative to not depend on initial energy, to allow nonzero initial_vy and initial_v_vertical:
    which_belong_here = false(length(trajectories_ungrouped.t),1);
    if mod(length(trajectories_ungrouped.t), source.ion_count) == 0 && length(trajectories_ungrouped.t) == source.ion_count * (length(run_index)-returned_runs+1)
      % OK, the same number of particles for all runs. (Assumes all runs (e.g. electrons and ions) were with same number of particles!)
      which_belong_here( 1:source.ion_count ) = true; % Because entries are removed from the trajectories_ungrouped when used, the first block of indices is always the one to select
      which_belong_here = which_belong_here & trajectories_ungrouped.charge == source.charge ...
                                            & abs(trajectories_ungrouped.mass - source.mass) < 5E-5;
    else
%       if length(which_belong_here) > 0
%         % Trajectories were given, but none could be selected. This seems to be an error.
%         warning('readsim3:trajectory_selection', 'Trajectories were given, but none could be selected. This seems to be an error...');
%       end
    end

    if sum(which_belong_here) ~= source.ion_count
      warning('readsim3:trajectories', sprintf('Failed to find the correct number of trajetories for summary %d in run %d: Got %d instead of %d. ', index, run_index(index), sum(which_belong_here), source.ion_count)); % NOTE: this may repeat message to concatenate for several cases
      trajectories = struct('y',[], 't', [], 'mass',[], 'charge',[], ...
        'energies',[], 'elevation',[], 'z0', [], 'y0', [], 'z', [], 'x',[], 'x0',[], 'azimuth',[] );
    else
      trajectories = struct('y',trajectories_ungrouped.y(which_belong_here), 't', trajectories_ungrouped.t(which_belong_here), ...
        'mass',trajectories_ungrouped.mass(which_belong_here), 'charge',trajectories_ungrouped.charge(which_belong_here), ...
        'energies',trajectories_ungrouped.energies(which_belong_here), 'elevation',trajectories_ungrouped.elevation(which_belong_here), ...
        'z0', trajectories_ungrouped.z0(which_belong_here), 'y0', trajectories_ungrouped.y0(which_belong_here), ...
        'z', trajectories_ungrouped.z(which_belong_here), 'x',trajectories_ungrouped.x(which_belong_here), ...
        'x0',trajectories_ungrouped.x0(which_belong_here), 'azimuth',trajectories_ungrouped.azimuth(which_belong_here) );
      % Remove the used particles from the ungrouped set
      trajectories_ungrouped.y(which_belong_here) = [];
      trajectories_ungrouped.t(which_belong_here) = [];
      trajectories_ungrouped.mass(which_belong_here) = [];
      trajectories_ungrouped.charge(which_belong_here) = [];
      trajectories_ungrouped.energies(which_belong_here) = [];
      trajectories_ungrouped.elevation(which_belong_here) = [];
      trajectories_ungrouped.z0(which_belong_here) = [];
      trajectories_ungrouped.y0(which_belong_here) = [];
      trajectories_ungrouped.z(which_belong_here) = [];
      trajectories_ungrouped.x(which_belong_here) = [];
      trajectories_ungrouped.x0(which_belong_here) = [];
      trajectories_ungrouped.azimuth(which_belong_here) = [];
    end
  end
  
  
  % Save it all in one structure
  fly = struct('source',source, 'rays',rays, 'ray_columns', [], ...
    'misses', source.ion_count-sum(rays(:,2)), 'overall',overall, ...
    'trajectories',trajectories, ...     % no longer: 'recorded_data',ion_recording, 
    'electrode_index', electrodes, 'potentials',potentials, 'V', V, ...
    'lens_parameters_A',lens_parameters_A, 'lens_parameters_B',lens_parameters_B, 'imaging_dimensions',imaging_dimensions, ...
    'geometry', geometry);
    % Since geometry info is only available with some new Lua scripts, users should always 
    % first verify that  isfield(fly, 'geometry') && ~isempty(fly.geometry)  returns true.
  
  fly.ray_columns = columns; % Can't assign cell directly in struct-consructor, it would make multi-element struct.

	% Update the potential energy variables (and tweak ratio) so they don't
  % have opposite sign when lens is made for negative particle.
  fly.lens_parameters_A.U_ext   = fly.lens_parameters_A.U_ext   * charge;
  fly.lens_parameters_A.U_bend  = fly.lens_parameters_A.U_bend  * charge;
  fly.lens_parameters_A.U_tweak = fly.lens_parameters_A.U_tweak * charge;
  fly.lens_parameters_B.U_ext   = fly.lens_parameters_B.U_ext   * charge;
  fly.lens_parameters_B.U_bend  = fly.lens_parameters_B.U_bend  * charge;
  fly.lens_parameters_B.U_tweak = fly.lens_parameters_B.U_tweak * charge;
  %(no longer needed, defined via U_ext) fly.lens_parameters_A.tweak1  = fly.lens_parameters_A.tweak1  * charge;
  %(no longer needed, defined via U_ext) fly.lens_parameters_B.tweak1  = fly.lens_parameters_B.tweak1  * charge;

  
  if returned_runs == 1
    flys = fly; % Create the "flys" variable (for struct-array if more runs will be returned)
  else % Append to existing list of runs
    flys(returned_runs,1) = fly;
  end

end
if returned_runs == 0
%   error('No runs were selected in %s.', filename)
  warning('No runs were selected in %s.', filename)

%   % Define an empty output anyway (previously it seems NaN was returned, but caller not prepared to handle it)
%   flys = struct('source',{}, 'rays',{}, 'ray_columns', {}, ...
%     'misses', {}, 'overall',{}, 'trajectories',{}, 'electrode_index', {}, 'potentials',{}, 'V', {}, ...
%     'lens_parameters_A',{}, 'lens_parameters_B',{}, 'imaging_dimensions',{} )
end
clear s file_contents fly

if ~isempty(flys(1).trajectories.t) || ~isempty(flys(end).trajectories.t)
  % Now that energy filtering is not use when picking the trajectories, make a check that the ordering seems correct.
  % Although randomization is active, it should be re-seeded and reproduced for each round (energy), so I expect monotonicity to hold.
  energy_hint = NaN(length(flys), 5); %[mass,nominal_without_initial_velocity_offset, total_min,total_max,total_mean]
  for i = 1:length(flys)
    energy_hint(i,:) = [flys(i).source.mass, flys(i).source.energy, minmax(flys(i).trajectories.energies'), mean(flys(i).trajectories.energies)];
  end
  for m = unique(energy_hint(:,1))'
    % For each mass (Since a constant velocity corresponds to different kinetic energies for different masses,
    % and mass-check already done in trajectory selection, only order within mass needs to be checked here)
    energy_hint_for_m = sortrows(energy_hint(energy_hint(:,1)==m,:), [2]);
    if sum(any(diff(energy_hint_for_m,[],1) < -0.04)) > 0
      % If out-of-order (by more than 0.04eV) among the energy_hint appears in any column, then warn
      warning('readsim3:trajectory_selection', 'After loading %d runs with, the energies within them appear to be misordered for mass %.3f u. This could indicate a failure of the de-mixing procedure for trajectories_ungrouped. %d of 3 columns disagreeing, if <3 then maybe it is an artefact when using initial_vy or initial_v_vertical.', length(flys), m, sum(any(diff(energy_hint_for_m,[],1) < -0.04)));
    end
  end
end

% Inner helper function to parse a line that is repeated for side A and side B
function [lens_parameters, potential_string] = read_lens_parameters(potential_string, lens_parameters)

[values,count,~,nextindex] = sscanf(potential_string, ['%*[^:]: E_ext=%g V/cm, E_tube=%g V/cm, %*[bend]=%g, t%*[^=]=%g tw''=%g, ext/rep%g'], 6); % try reading "tw'", but accept also if not succeding
if count >= 4
  if ~isfield(lens_parameters, 'bend') || isempty(lens_parameters.bend) || isnan(lens_parameters.bend)
    lens_parameters.bend = values(3);
  else
    % keep value calculated from potentials, since it keeps more significant digits than the value printed in log
  end
  if ~isfield(lens_parameters, 'tweak1') || isempty(lens_parameters.tweak1) || isnan(lens_parameters.tweak1)
    lens_parameters.tweak1 = values(4); % = V_Atweak/W_ext
  else
    % keep value calculated from potentials, since it keeps more significant digits than the value printed in log
  end
  

  % NOTE: old V_TUBE and V_TWEAK names were changed to U_bend and U_tweak for consistency, because
  % their sign is such that they represent potential energy difference = (electric voltage) * charge.
  if count >= 5
    lens_parameters.tweak2 = values(5); % The variant 2 of tweak parameter: voltage with respect to untweaked voltage (considering plateau, if any).
    lens_parameters.tweak3 = NaN;%values(6); % Skipping tweak3 (based on interpolation, using Z-coordinate for lens electrode but unaware of plateau and no exact meaning)
  else
    lens_parameters.tweak2 = NaN;
    lens_parameters.tweak3 = NaN;
    nextindex = nextindex - 1; % Go back before the character that caused incomplete match (if count == 4 with the tweak2= pattern)
  end
  if count >= 5
    lens_parameters.extractor_repeller_ratio = values(6);% new for VMI-style lens comparison to litterature
    % This is essentially (for a two-field spectrometer, with last electrode at ground)
    % a different parameterization of the bend ratio:
    %   bend = 2 * extractor_ratio / (1 - extractor_ratio)
    %   extractor_ratio = bend/(bend + 2)
  else
    lens_parameters.extractor_repeller_ratio = NaN;
  end
  lens_parameters.E_extraction = values(1);
  lens_parameters.E_tube = values(2);

  nextindex = nextindex-1 + strfind(potential_string(nextindex:min([length(potential_string) nextindex+200])), sprintf('\n')); % skip past the used info, then to the end of this line
  potential_string = potential_string(nextindex(1):end);

else
  if ~isfield(lens_parameters, 'bend')
    lens_parameters.bend = NaN;
  end
  if ~isfield(lens_parameters, 'free')
    lens_parameters.free = NaN;
  end
  if ~isfield(lens_parameters, 'plateau')
    lens_parameters.plateau = NaN;
  end
  if ~isfield(lens_parameters, 'tweak1')
    lens_parameters.tweak1 = NaN;
  end
  lens_parameters.tweak2 = NaN;
  lens_parameters.tweak3 = NaN;
  lens_parameters.extractor_repeller_ratio = NaN;
  lens_parameters.E_extraction = NaN;
  lens_parameters.E_tube = NaN;
end
