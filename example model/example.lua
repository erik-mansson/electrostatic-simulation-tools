-- Script for a SIMION workbench (.iob & .pa0) by Erik P. Månsson.
-- Lund University 2010-2014. IFN-CNR/Politecnico di Milano 2014-2017. FS-ATTO CFEL DESY 2017-2020.
-- 
-- This file should be named like the model (.pa0) file but with the .lua extension, and some constants
-- need to be defined suitable for the model's geometry, down to the source_point_y, source_point_z, 
-- VERTICAL_CENTRE and Y_CENTRE coordinates, while the bulk of the code and most variables are model-independent.
-- 
-- Since the program has a long history, initially with featuers for usage interactively
-- in the SIMION GUI but later more used in command-line calls by runsim3.m and readsim3.m,
-- there are many variables that are ignored with the default settings. Limitations in the
-- availble invocation points from SIMION (e.g. how to know when a particle is the last particle 
-- of its energy?), and starting to develop this in earlier SIMION versions, have also contributed 
-- to a complex structure of the program's statistics-gathering and summary-printing parts.
-- 
-- SIMION (https://simion.com) is used to simulate electron and/or ion trajectories in electrostatic
-- spectrometers. This script was developed for designing and tuning spectrometers intended for 
-- 3D momentum imaging (time of flight + 2D position), velocity map imaging (VMI in two dimensions) 
-- as well as mass resolution (but the script can not directly show a mass resolution, only raw data
-- and statistics to be judged manually or by a higher-level program).
-- 
-- Key functionalities:
-- 
-- * The particles' initial position can be set in a deterministic pattern (e.g. line or grid-like)
--   specified by integer values (>=0) for source_point_pattern. A negative value (which is the default)
--   instead leaves the position defined by SIMION GUI under Particles -> Define particles, which in turn
--   is loaded from the .fly2 file if there is one (e.g. a gaussian3d_distribution is used in the example.fly2).
--   Search for source_point_pattern in this script to learn about the specific settings and further
--   variables used, e.g. _source_point_spacing, source_point_y and source_point_z.
-- 
-- * The particles' initial direction (but not energy) is overridden based on the _isotropic variable,
--   and this also controls how the directions are grouped for in a summary printing the average
--   coordinates, the standard deviation (small for good VMI-focus, large if initial position effect remains),
--   and summarized statistics in the form of percentages to easily see in GUI whether a change makes
--   the VMI-focus better or worse. The most used settings are _isotropic=2 for 2D-VMI and =0 for 3D-momentum imaging.
--   To use random directions but group them all without caring about direction (e.g. to just judge mass resolution)
--   _isotropic=1 together with _number_of_directions=1 can be useful too.
--   Either way, the kinetic energy of particles is defined in the GUI (initialized from loaded .fly2-file).
-- 
-- * The potentials of the electrodes (the voltages) can be defined in different ways:
--   DoNotSetPotentials=1 (default) leves the potentials as defined in SMION GUI (PAs -> Fasst Adjust Voltages...).
--   DoNotSetPotentials=0 uses many variables in the script to update the potentials, e.g. including two resistor 
--     chains and a few "free" potentials. (e.g. _U_Afirst, _U_Alast, _U_AMCP, _U_Afree, _A_free, A_FIRST and A_LAST)
--   DoNotSetPotentials=-1 uses a given array with all the potentials, in the global variable _G.U (or the variable U). 
--     This is very useful for models with more complex electrode arrangement than a pair of simple resistor chains.
--     When tha MATLAB function runsim3.m launches a simulation by the command line, it writes a custom .fly2 file
--     defining particles as well as _G.U. For use in the GUI, you can adjust all voltages onelines like the following
--     on the SIMION command line, e.g. for a model with three potentials one could type:
--      _G.U={-4000,-4000*0.7,-50}; _G.U[0]=0; U=_G.U; adjustable DoNotSetPotentials=-1;
-- 
-- * For the printing of statistics-summaries to work (correct grouping of directions and energies),
--   consistent values are necessary in several variables. Basically, the global variable _G.FLY2_IONS_PER_GROUP 
--   should be equal to the "Num particles" in the GUI's Particles -> Define particles, and this value
--   should be an integer multiple of the adjustable variable _number_of_directions. Each direction gets 
--   _G.FLY2_IONS_PER_GROUP/_number_of_directions particles, e.g. 32/8=4 by default. If you change the 
--   GUI's "Num particles" (in particle definition GUI) you should use the command line to set _G.FLY2_IONS_PER_GROUP
--   to the same value. (Making "Num particles" a multiple of _G.FLY2_IONS_PER_GROUP may work to print multiple 
--   summaries of particles with the same energy, instead of one big group with better statistics.)
--
-- * SIMION can be configured to log details about every single particle trajectory, to allow re-analyzing
--   the resolutions in more advanced ways (including mass resolution) in another program. When called by runsim3.m, 
--   The options --recording=record3.rec --recording-output=out.txt means that the format specifier record3.rec
--   specifies the logging format and that output is written to out.txt (including raw trajectory log and 
--   this script's summary printouts)
--
-- * The individual particle impact times and positions can also be written to a file in the DLT file
--   format also used for experimental 3D-momentum data recorded via a delay-line detector (DLD).
--
-- * If a file called adjust.lua exists, it will be run too at startup, as another way of overriding many variables.
--   Since the use of a fixed filename isn't  suitable for multi-threaded simulations
--
-- NOTE: In SIMION the cylinder axis is named x, while z is preferred for this symmetry axis in non-SIMION 
-- descriptions of spectrometers. Most of the script's variables thus use x and z swapped compared to the
-- actual SIMION coordinate system.
-- By convention, y is used for the optical axis, along which the initial points (ionization volume)
-- has the largest length (assuming laser if focused to a narrow beam waist).
-- Due to cylidrical symmetry (at least of VMI-spectrometers), y can also be thought to represent the radius
-- (orthogonal to cylinder axis), while in the general model y is only one of the components of r.
-- The coordinate that is neither spectrometer symemtry axis nor optical axis (i.e. Z in SIMION, x in variable naming)
-- is not used for velocity directions in the typical settings for _isotropic, only for (Gaussian)
-- size of the ionizatiaon volume (spread in source points) and for an optional initial_v_vertical 
-- velocity to represent a common jet-velocity offset, of relevance mainly for ions with longer times of flight. 
-- I.e. a particle with kinetic energy in GUI defined as zero gets initial_v_vertical as vertical velocity,
-- and particles defined with nonzero kinetic energy gets the combined velocity vector.
--
simion.workbench_program()

-- 2013-11-04 setting trajectory quality to something higher than default (3).
-- (9 is max for automatic steps, but 100-109 range can be used for always 0.1mm max, http://simion.com/info/faq.html#how-does-the-trajectory-quality-factor-t-qual-relate-to-time-step-sizes)
sim_trajectory_quality = 6 -- maybe this was a bit slow for running thousands of simulations? Restored 2015-09-11 18:15
--sim_trajectory_quality = 4 -- used for optimizations until 2015-09-11

-- These _U... and _A... parameters are only used if DoNotSetPotentials is set to 0. Useful for simple models and adjusting voltages from GUI.
adjustable _U_Bfirst  = -4000 -- The potential [V] applied to electrode 1 (most distant from detector)
adjustable _U_Afirst  = -2800 -- The potential [V] applied to electrode 8 (the separating mesh). (Ignored if R_control_8==0)
adjustable _U_Bfree = -10 -- The potential [V] applied to electrode _free_index (main Einzel lens electrode on A-side, isolated from other electrodes)
adjustable _U_Afree = -0 -- The potential [V] applied to electrode _free_index (main Einzel lens electrode on A-side, isolated from other electrodes)
adjustable _U_Blast = -10 -- The potential [V] applied to electrode B_LAST (closest to phosphorescent screen).
adjustable _U_Alast = 0 -- The potential [V] applied to electrode A_LAST (closest to DLD).
adjustable _A_free = 1 -- which electrode is isolated from the others on A-side, the main Enzel lens electrode, given as the offset from A_FIRST
adjustable _A_plateau = 0 --how long is the plateau (short circuited electrodes) on the A-side, given as offset from A_FIRST
adjustable _U_BMCP = -10   -- Potential [V] at the "B"-MCP front 
adjustable _U_AMCP = 0   -- Potential [V] at the "A"-MCP front 

-- If the design is modified, these parameters make it simple to adapt the program.
-- The electrode indexing is described by the following constants, dependent on the SIMION potential array:
--   (Values in parentheses give the index used in the "doublesided_20110117" workbench.)
--   WALLS (1) is the chamber walls, always kept grounded. Can be set to 0 if the workbench has no walls.
--   B_MCP (2) is the B-MCP (B is the phosphorescent screen side, "front" is towards source point)
--   B_LAST (6) is the electrode closest to the B-MCP, with a meash in this case.
--   B_FIRST (8) is the extractor electrode towards the B-side.
--     The sequence of indices from B_LAST to B_FIRST (6,7,8) belong to the B-side (while 3,4,5 were nonexistent in doublesided_20110117).
--     The sequence of indices B_FIRST-B_FREEGROUP_COUNT to B_FIRST-1 (7) act as an Einzel-like tweak-electrode on the A-side, with a common potential _U_Bfree.
--   From B_FIRST to A_FIRST (8,9,10,11) are the electrodes in and meshes surrounding the source region.
--   A_FIRST (11) is the extractor electrode towards the A-side, starting a sequence in indices up to A_LAST.
--   A_LAST  (25) is the electrode closest to the B-MCP, with a meash in this case.
--     The sequence of indices A_FIRST to A_LAST (12,...,25) are the electrodes on the A-side.
--     The sequence of indices A_FIRST+1 to A_FIRST+_A_free (adjustable variable!) act as an Einzel-like tweak-electrode on the A-side, with a common potential _U_Afree.
--     The sequence of indices A_FIRST+1 to A_FIRST+_A_plateau (adjustable variable!) except those in the Einzel-electrode, are kept at the same potential as A_FIRST.
--   A_MCP (26) is the A-MCP front  (A is the delay-line detector side, "front" is towards source point)
WALLS = 0; --# TODO make walls
B_MCP = 1; B_LAST = 1; B_FIRST = 1;
A_FIRST = 2; A_LAST = 3; A_MCP = 3; --(currently no A-MCP, but 0 not allowed)
B_FREEGROUP_COUNT = 0; -- this is not an index, but the number of plates that form a free electrode group on the B-side


-- The ion start positions (except along spectrometer symmetry axis) are set by the program.
adjustable _number_of_directions = 8 -- The number of directions must be a factor in the number of ions flown (.fly2-file or GUI setting). For _isotropic=0: typically 176 ions for 16 directions, or 88 for 8; For _isotropic=2: typical number of direction are 16, 8 (even numbers avoid theta=0), 9, 5 or 3.
adjustable _isotropic = 2 -- 0: discrete directions 0 to 360 degrees; 1: isotropic, random distribution, 2: 2D-VMI-style directions +90 to -90 degrees, 3&4: VMI-style with approximately (discretized, deterministic) isotropic 3D distribution.
adjustable _source_point_spacing = 0.8 -- was 1 with pattern=0. -- Was 0.5 even earlier. --[mm] (changed from 1 2009-09-24, changed to 1 2010-10-08)
adjustable source_point_pattern = -2 -- <0: use SIMION FLY or FLY2 configuration, 0: elongated only along r|y (standard before 2010-10-14), 1: both z|x and r|y, 2: both in better patten (assymetric with respect to y-sign), ... (see code)
-- The source point is for source_point_pattern==-1 given by .FLY2 file (command line usage) or user interface's Particle->Define (defaulting to a fly2-file for the workbench), 
-- but for non-negative source_point_pattern the following coordinate is used:
adjustable source_point_y = 0 --[mm] Centre of distribution of particle starts (if source_point_pattern>=0), along optical axis. See also detectorA_y and detectorB_linecoeff[4]. Old name: source_point_offset_y (same meaning)
adjustable source_point_z = 90 --[mm] - Centre of distribution of particle starts (if source_point_pattern>=0), along spectrometer axis (x in SIMION).

adjustable initial_vy = 0 --[mm/us]=[km/s] Initial offset velocity
adjustable initial_v_vertical = 0 --[mm/us]=[km/s] Initial offset velocity, to simulate skimmed pulsed supersonic jet. SIMION z axis (NOTE: this was previously called initial_vz which was a confusing use of SIMION-coordiates instead of the swapped x/z as for most variables)

-- Allows you to keep any potentials set using "fast adjust" in Simion, while running the program
adjustable DoNotSetPotentials = 1
_G.ions_flown = 0 -- To find out how many ions were flown (per group) (not to be adjusted by user)
_G.ions_flown_all_groups = 0 -- To find out how many ions were flown globally (not to be adjusted by user)

-- The number of source points is determined by
-- <number of ions to be flown (user/.FLY-specificed)> / _number_of_directions
-- Or the inverse:
-- <number of ions to be flown> = _number_of_directions * (1 + <desired_y_range>/_source_point_spacing)
-- E.g. 8 directions, from y=-5 to y=5 spaced by 1 [mm] gives <number of ions> = 88
adjustable detectorA_x = 283 --[mm] A-MCP front surface position (depends on the geometry). This is actually "Z" coordinate in typical coordinate system (along spectrometer axis), but x in SIMION.
detectorA_y = 0 --[mm] Central y-coordinate of A-MCP. Introduced 2015-12-30 to allow proper VMI analysis in Lua also in the new designs where detectors are not at y=0. See also source_point_y.
adjustable detectorB_x = 49  --[mm] B-MCP front surface centre position (depends on the geometry). This is actually "Z" coordinate in typical coordinate system (along spectrometer axis), but x in SIMION.
--(both4) detectorB_linecoeff = {-5.1446, -1, 2059.6177, 187} -- {k, -1, -d, MCP_centre_y}: k=-1/tan(det.angle), d = k*MCP_centre_x - 1*MCP_centre_y; using MCP imported as separate PA
detectorB_linecoeff = {-0, -0, -999, -999} -- {k, -1, -d, MCP_centre_y}: k=-1/tan(det.angle), d = k*MCP_centre_x - 1*MCP_centre_y; using MCP imported as separate PA

Z_EXTRACTOR_A = 39.75 --[mm] coordinate along spectrometer axis for the first mesh on A-side
Z_EXTRACTOR_B = 24.75 --[mm] coordinate along spectrometer axis for the first mesh on B-side
Z_END_A = 999 --[mm] coordinate along spectrometer axis for the last mesh on A-side
Z_END_B = 1.5 --[mm] coordinate along spectrometer axis for the last mesh on B-side
Z_FREE_B = 54.75 --[mm] rough coordinate along spectrometer axis for the free electrode on B-side (to define a default value when not much lensing)
DZ_PLATE_THICKNESS = 1 --[mm] thickness of a plate
DZ_PER_INDEX_A = 10 --[mm] distance between plates + plate thickness (if variable, should match the A-side between Z_EXTRACTOR_A and the free electrode)
R_GEOMETRY_MAX = 56 --[mm] prevent hangup if source points would be outside simulation volume
Z_GEOMETRY_MAX = 400 --[mm] (along TOF axis) attempt to prevent missing printout for particles leaving simulation volume (if no outer ground-electrode stops them)
VERTICAL_CENTRE = 0 --[mm] for the non-vertically symmetric "both8eu" but 0 in all other/earlier models
Y_CENTRE = detectorA_y --[mm] so far, detector A has always been centered so we can reuse that value

--------- THE PART FROM HERE AND BELOW IS THE SAME FOR ALL WORKBENCHES (until the simulation volume size is changed)

adjustable multipole_order = 1 -- Legendre Polynomial order for _isotropic==1. 1: dipole, 3: third order, 4: fourth order.
output_filename = "sim.dlt" -- where a DLT file will be written (if enabled)
output_DLT_file = nil -- will hold the file handle
adjustable write_DLT = 0 -- 0: no, 1: head and hits, 2: hits only, 3: hits and foot, 4: all
adjustable time_offset_DLT = 0 -- [us] (microseconds, the SIMION time unit) offset for the TOF
if _G.DLT_comment == nil then --when not running as invoked from commandline.lua
  _G.DLT_comment = ""
  _G.DLT_hit_count = 0
  _G.DLT_trigger_count = 0
end
adjustable random_seed = 20150319 -- To make pseudorandom particle distributions reproducable.

-- Since sim_update_pe_surface only seems possible to set from other_actions
adjustable potential_display_needs_updating = 1  --(only used by the program, not for user input)
adjustable has_shown_potentials = 0 --(only used by the program, not for user input) to do it only once

-- Re-define some adjustables to other values than in the main/default program.
CMD,errmsg = loadfile('adjust.lua');
if CMD then
  setfenv(CMD,getfenv()); CMD()
  print('Adjustments loaded. Press the shortcut keys: V, F, C and A.\n  (Tip: maximize the window and set the value 50 at the Auto before pressing A)')
else
  print('No external adjustments.')
end

math.randomseed(random_seed) -- To make pseudorandom distributions generated in Lua reproducable (not sure if this is re-run for each fly or only once, but needed when called from command line/Matlab)
 
-- To fly multiple energies and print separate statistics for them,
-- a FLY2 file can be prepared to set the global variable G.FLY2_IONS_PER_GROUP
-- (SIMION 8.0 runs all particles together, not finishing first group before starting second, so effort is required to not mix the data)
-- NOTE: SIMION 8.1 has better ways to control iteration by function segment.flym()...
if _G.FLY2_IONS_PER_GROUP == nil then
  _G.FLY2_IONS_PER_GROUP = 2147483648 -- put a large number instead (more ions than 2^31 won't be used)  
end
  
 -- When computing relative resolutions, the standard deviations of each bunch (direction) is used. However
 -- for very narrow bunches a lower limit that represents the experimental resolution (detector bin size and random fluctuations during measurement).
 -- These thresholds are named RES_... here. 
-- Until 2013-11-04 the absolute resolutions were RES_R_ABS = 0.5, RES_R_REL = 1.0, RES_TOF_ABS = 0.8e-3, RES_TOF_REL = 1.0e-3
-- Seems they can be lowered now (new values are estimated FWHM but maybe std assumed in this program, still a reduction from old).

RES_R_ABS  = 0.25 -- [mm] For averaging of absolute radial resolution ("tot sy").
RES_R_REL   = 0.5 -- [mm] For averaging of relative radial resolution ("<sy/y>" and "<1/R>"). 1 ns in DLD corresponds to roughly 0.5 mm.
                  -- Maybe a bit small but since the bunches have roughly rectangular probability distribution and not gaussian the std. is "small" compared to the real data,
                  -- so I think 0.5 mm resolution referes more to gaussian and that a smaller value could be possible
                  -- (and small value is in line with version before 2010-09-10 which did not limit the radial std that was used to get relative radial resolution).
                  -- New: Actually, DLD wires are probably 1 mm apart, so raising it again.
MIN_R_REL  = 1    -- [mm] For computing relative radial resolution, don't divide by radius^2 smaller than this (relative error is undefined if y=0).
                  -- The (relative error)^2 is weighted by by MIN_R_WEIGHT + abs(sin(elevation) in the RMS so that
                  -- the error is (almost) ignored at elevations 0 and 180 degrees.
MIN_R_WEIGHT = 0.001 -- The (very small) weight used at 0 and 180 degrees (see above),
                  -- where the "relative error" (after thresholds) often equals RES_R_REL/MIN_R_REL (=40% with current values)

HIGH_R_REL = 30  -- [mm] For computing relative radial resolution for directions (elevations) near 0 and 180 degrees.
                     -- Since the radial spread does not blur energy info here (only angular info or overlap with other features)
                     -- it is "almost ignored" by dividing the std. (or RES_R_REL) with a "large radius" HIGH_R_REL.
                     -- (current values gives at least 0.4/30=1.3% for the 0 and 180 degree directions.)

RES_TOF_ABS = 0.3e-3 -- [us] For averaging of absolute TOF resolution ("st")
                     -- RES_TOF_ABS doesn't affect the "overall" "<1/R>" resolution presented.
RES_TOF_REL = 0.5e-3 -- [us] For averaging of relative TOF resolution (forward-backward direction difference, "<st/t">).
                     -- Maybe TOF detector can do a bit less than 1 ns, but quality of data would be worse and "too high" acceleration
                     -- would be favored too much (since also radial resolution often gains from it), by having a somewhat high
                     -- threshold in RES_TOF_REL this hopefully balances it a bit. Also the simulated "bunch std." is not a good model
                     -- of experimental source volume and inhomogenities, so attempting to minimize the simulated "buhch std." near detector
                     -- resolution doesn't seem relevant.


-- Variables used to calculate beam widths at detector.
-- Indexed by [group_index][direction_index], starting at 1.
local det_count = {}  --arrays in Lua are simply numerically indexed tables, and grow as needed
local det_y_min = {} 
local det_y_max = {}
local det_y_mean = {}
local det_y_Q = {} -- for computation of standard deviation
local det_y2_mean = {}
local det_y2_Q = {} -- for computation of standard deviation
local det_t_mean = {}
local det_t_Q = {}-- for computation of standard deviation
local copy_of_ion_number = -1
local highest_group_index = -1
_G.ion_energy = {} -- Ion initial energy, indexed by group_index, starting from 1
local group_charge = {} -- charge by group_index, starting from 1
local group_mass = {} -- charge by group_index, starting from 1


------- Additional library/helper functions:
-- "iif(c, a, b)" is a replacement for C-style "c ? a : b":
function iif(condition, when_true, when_false)
  return (condition and when_true) or when_false
end
-- Round number to string with given precision
function rounds(number, decimal_precision)
  return string.format("%." .. (decimal_precision or 0) .. "f", number)
end
-- An improved alias for print(string.format( ))
function printf(...)
  -- Improve the string.format so that zeros in exponential format are removed, e.g. "2.2e+7" instead of "2.2e+007"
  local s = string.gsub(string.format(...), "(%d)e([-+])(0+)(%d)", "%1e%2%4")
  return print(s)
end
-- An alias for error(string.format(...))
function errorf(...)
  local s = string.gsub(string.format(...), "(%d)e([-+])(0+)(%d)", "%1e%2%4")
  return error(s, 2)
end
function math.max_abs_tbl(t)
  local m = nil
  for _,v in pairs(t) do
    if m then
      m = math.max(m, math.abs(v))
    else
      m = math.abs(v)
    end
  end
  return m
end


TWO_PI = 6.283185307179586476925286766559;
DEGREE = 3.1415926535897932384626433832795/180;
-- For DLT-file writing:
HW_TICK = 0.000025 --[us] hardware clock period time = 25ps
-- The factors that convert between distance on detector and
-- delay-line-time coordinate [nm/ns]. r = G * t <==> t = r / G
G_x = 2 * 550.0 --2 * G_x [mm/us] is SIMION units. (was 550000.0 [m/s]=[um/us] in LabVIEW)
G_y = 2 * 520.0 --2 * G_y [mm/us]. The factor 2 is included to avoid separate +-t_x/2 and +-t_y/2 in the file writing expressions
-- Convert a 24-bit (unsigned or signed) number (integer) to a 3-byte string.
-- NOTE: there is no overflow/underflow checking, values simply wrap.
function u24_as_string(value)
  if value < 0 then
    --error "Negative value can't be handled by u24_as_string."
    value = 16777216 + value -- 2^24 + value gives the expected string (if -2^23<=value<0)
  end
  local t = math.floor(value / 256);
  return string.char( math.mod(math.floor(t / 256), 256), math.mod(t, 256), math.mod(value, 256));
end
-------

-- Temporary array used during calculation of potential values, for storage until assignment
-- into adj_elect[] by init_p_values (Don't know why this was introduced, can't they be applied directly?)
-- U is indexed just as adj_elect (unlike in previous Lua programs).
if DoNotSetPotentials ~= -1 then
  U = {}
end

-- Calculate the electic potential values, store in U[], but not yet in adj_elect[]
function segment.initialize()
  if has_shown_potentials == 0 then -- only first time
    local U_Afree_if_nonlensing_v2 
    if DoNotSetPotentials == 0 then
      -- Using adjustable variables to control potentials.
      -- NOTE: the code in this block depends on how the circuit (resistors) are connected.
      -- Currently A-side is based on old "simple7b" with adaptations.
      -- B-side has the first two electrodes acting as a single free (at same potential) and also as the plateu,
      -- so the voltage steps are equal if this free plate-pair is ignored.

      U[WALLS] = 0; -- always grounded chamber
      U[A_MCP] = _U_AMCP;
      U[B_MCP] = _U_BMCP;

      local A_COUNT = A_LAST - A_FIRST + 1; -- for A the indices increase from FIRST to LAST
      local B_COUNT = B_FIRST - B_LAST + 1; -- for b the indices increase from LAST to FIRST

      ---- Conditions for the A-side: ----
      if _A_plateau < 0 or _A_plateau > A_COUNT-1 then
        errorf("_plateau_to_index %g is out of the allowed range from 0 to %d.", _A_plateau, A_COUNT-1)
      end
      if _A_free < 0 or _A_free > A_COUNT-1 then --(Allowed to use =0 if no free potential)
        errorf("_A_free %g is out of the allowed range from 0 to %d. (0 may be used to not use any free potential.)", _free_index, A_COUNT-1)
      end
      -- There are still two resistors between _free_index-1 and _free_index+1 but
      -- electrode _free_index is isolated (free) from the point where the resistors are connected.
      -- Thus, all the potentials may be calculated as if there was no free electrode, and then the free electrode
      -- potential is simply changed (works also when _free_index < _plateau_to_index).
      local number_of_drift_resistors = A_COUNT - 1 - _A_plateau
      local i
      
      ---- Determine potentials in middle region ("acceleration region") ----
      for i = B_FIRST, A_FIRST do
        U[i] = _U_Afirst + (_U_Bfirst-_U_Afirst) * (A_FIRST-i)/(A_FIRST-B_FIRST)
      end
      if math.abs(U[B_FIRST] - _U_Bfirst) > 1e-5 then
        errorf("Precision error in U_Bfirst result: %19g instead of %19f.", U[B_FIRST], _U_Bfirst)
      end

      ---- A-side (to DLD)  ----
      -- Determine potentials in plateau
      for i = A_FIRST+1, A_FIRST+_A_plateau do
        U[i] = _U_Afirst -- short circuit to electrode A_FIRST
      end
      -- Determine potentials after plateau, via the current that flows trough this part of the circuit
      for i = 1, number_of_drift_resistors do
        U[A_FIRST+_A_plateau+i] = _U_Afirst + (_U_Alast - _U_Afirst) * i / number_of_drift_resistors
      end
      U_Afree_if_nonlensing_v2 = U[A_FIRST+_A_free]
      if _A_free > 0 then
        -- Finally, override the potential for the free (lens) electrode that is isolated from "its normal position" in the resistor chain
        -- Changed circuit 2011-09-03: Free electrode is everything from first after mesh to the index specified (2 in circuit now)
        for i = 1, _A_free do
          U[A_FIRST+i] = _U_Afree -- isolated from others and externally controlled
        end
      end
      if math.abs(U[A_LAST] - _U_Alast) > 1e-5 then
        errorf("Precision error in U_Alast result: %19g instead of %19f.", U[A_LAST], _U_Alast)
      end

      ---- B-side (to phosphorescent screen) ----
      for i = B_FIRST-B_FREEGROUP_COUNT, B_FIRST-1 do --NOTE: this loop typically runs zero or one iterations (with B_FREEGROUP_COUNT being 0 or 1)
        U[i] = _U_Bfree
      end
      -- Using the same resistance between all non-free electrodes
      number_of_drift_resistors = B_COUNT - 1 - B_FREEGROUP_COUNT
      for i = 1, number_of_drift_resistors do
        U[B_FIRST-B_FREEGROUP_COUNT-i] = _U_Bfirst + (_U_Blast - _U_Bfirst) * i / number_of_drift_resistors
      end
      if B_COUNT > 1 and math.abs(U[B_LAST] - _U_Blast) > 1e-5 then
        errorf("Precision error in U_Blast result: %19g instead of %19f.", U[B_LAST], _U_Blast)
      end
      
      -- There is no resistor chain on the B-side now, so the default potential is not defined.
      -- Use a linear rough interpolation:
      local U_Bfree_if_nonlensing_v3 = ((Z_FREE_B-Z_EXTRACTOR_B)*U[B_LAST] + (Z_END_B-Z_FREE_B)*U[B_FIRST])/(Z_END_B-Z_EXTRACTOR_B)
      -- Also do it for the A-side (for version 3 of the tweak parameter)
      Z_FREE_A = Z_EXTRACTOR_A + DZ_PER_INDEX_A * _A_free + DZ_PLATE_THICKNESS / 2 --[mm] position of the free electrode on the A-side
      local U_Afree_if_nonlensing_v3 = ((Z_FREE_A-Z_EXTRACTOR_A)*U[A_LAST] + (Z_END_A-Z_FREE_A)*U[A_FIRST])/(Z_END_A-Z_EXTRACTOR_A)
    
    elseif DoNotSetPotentials == -1 then
      -- For case when potentials set not by the adjustable-variables,
      -- but directly from the U-array (from adjust.lua script).
      if U == nil or #U == 0 then
        -- If assignment was done via .fly2 script rather than adjust.lua, then global _G.U is used rather than U.
        U = _G.U
        if U == nil or #U == 0 then
          error('DoNotSetPotentials was set to -1 but no array of potentials given in U or _G.U.');
        end
      end
        
          -- Prevent printout from showing false (non-affecting) parameter values, show NaN instead
      _U_AMCP = U[A_MCP];
      _U_BMCP = U[B_MCP];
      _A_free = 0.0/0.0;
      _A_plateau = 0.0/0.0;
      B_FREEGROUP_COUNT = 0.0/0.0;
    end

    has_shown_potentials = 1
    print(" ") -- empty line
    if DoNotSetPotentials == 0 or DoNotSetPotentials == -1 then
      -- If potentials set via GUI of U-array:
	  if #U < A_MCP or #U < A_LAST then
		error(string.format("Got too few (%d) voltages in U for DoNotSetPotentials==-1 for A_LAST=%d and A_MCP=%d.", #U, A_LAST, A_MCP)) 
	  end
      -- Print the potentials (NOTE: will skip the unused electrodes)
      local tmp = string.format("U_AMCP = %g, U_Alast = %g, U_BMCP = %g, U_Blast = %g, A_free=%g, A_plateau=%g", U[A_MCP], U[A_LAST], U[B_MCP], U[B_LAST], _A_free, _A_plateau)
      tmp = tmp:gsub("-1[.]#IND", "NaN") -- show NaN instead of -1.#IND (note %f and %g give different strings for NaN)
      print(tmp)
    
      for i = 1, #U do -- print all potentials in the order of the workbench's potential arrays
        if U[i] then -- skip some unused indices (there may be more workbench array entries than electrodes)
          if i < #U and U[i+1] then
            printf("U[%d] = %g,  V = %g", i, U[i], U[i] - U[i+1])
          else -- can't show any voltage to/from the last electrodes
            printf("U[%d] = %g,  V = NaN", i, U[i])
          end
        end
      end
      -- New in version 3 (put after potential list to allow some compatibility)
      tmp = string.format("Indices: A_MCP=%d, A_LAST=%d, A_FIRST=%d, B_MCP=%d, B_LAST=%d, B_FIRST=%d, B_FREEGROUP_COUNT=%g, WALLS=%d",
             A_MCP, A_LAST, A_FIRST, B_MCP, B_LAST, B_FIRST, B_FREEGROUP_COUNT, WALLS)
      tmp = tmp:gsub("-1[.]#IND", "NaN") -- show NaN instead of -1.#IND (note %f and %g give different strings for NaN) for B_FREEGROUP_COUNT when NaN
      print(tmp)
    end
    
    if DoNotSetPotentials == 0 then
      -- Continuation of using adjustable variables to control potentials.

      -- Print "lens parameters" as the are defined now
      -- the z-coordinates should be in [mm] and the potentials in [V].
      function print_lensparam(side, z_repeller,z_extractor,z_end, U_rep,U_ext, U_free, U_free_if_nonlensing_v2, U_free_if_nonlensing_v3, U_end)
        -- Voltage across the interaction region (between the first meshes)
        local V_meshes = U_rep - U_ext -- [V] positive if ions towards this side, negative if electrons
        -- Voltage from the interaction|source point to the first mesh.
        --  V_ext = U_needle - U_ext (= V_meshes / 2 if the source is centered in extraction region)
        local V_ext = ((source_point_z-z_extractor)*U_rep + (z_repeller-source_point_z)*U_ext)/(z_repeller-z_extractor) - U_ext
        local V_tube = (U_ext - U_end)
        
        -- Accelerating field in the interaction region:
        local E_ext = V_meshes / math.abs(z_extractor-z_repeller) * 10 -- [V/cm] positive if ions towards this side, negative if electrons
        -- Drift-tube field if there was no plateau or free plate:
        local E_tube = V_tube / math.abs(z_end-z_extractor) * 10 -- [V/cm] positive if ions towards this side, negative if electrons
        
        printf("%s: E_ext=%6.4g V/cm, E_tube=%6.4g V/cm, b=%6.3f, tw=%6.3f tw'=%6.3f, ext/rep %5.3f", side,
          E_ext, E_tube, V_tube/V_ext, -- version 2 of the bend parameter (fits the formulas for TOF better), version 1 was: E_tube/E_ext
          (U_free - U_ext) / V_ext, -- similar to version 1 of the "tweak parameter" but using V_ext instead of V_meshes for generality (very simple, not good if lens would be far from extractor)
          (U_free - U_free_if_nonlensing_v2) / V_ext, -- version 2 of the "tweak parameter" definition (if free is within plateau it is just 2 * the version 1)
          (U_ext-U_end)/(U_rep-U_end) -- VMI style extractor/repeller ratio (with respect to detector&tube, ignoring that lens may deviate)
        )
          --(U_free - U_ext) / V_meshes, -- version 1 of the "tweak parameter" definition (simpler, but maybe more coupled with the bend parameter)
          --(U_free - U_free_if_nonlensing_v3) / V_ext, -- version 3 of the "tweak parameter" definition (ignoring the plateau, using interpolated V_free_if_nonlensing)
      end
      print_lensparam("A (right)", Z_EXTRACTOR_B,Z_EXTRACTOR_A,Z_END_A, U[B_FIRST],U[A_FIRST], _U_Afree,U_Afree_if_nonlensing_v2,U_Afree_if_nonlensing_v3, U[A_LAST])
      print_lensparam("B (left) ", Z_EXTRACTOR_A,Z_EXTRACTOR_B,Z_END_B, U[A_FIRST],U[B_FIRST], _U_Bfree,            U[B_FIRST]  ,U_Bfree_if_nonlensing_v3, U[B_LAST])
    end
  end

  if ion_number == 1 then -- Only once per fly (not again when using multiple energies or masses)
    if _isotropic == 1 then
      if _number_of_directions == -2 then
        -- Special case: simulate dipole transition where angular density is proportional to |y|^2 (Legendre polynomial of order 1, squared)
        print("# Using y-dipole distribution, random directions.")
      elseif _number_of_directions == 2 then
        -- Special case: simulate dipole transition where angular density is proportional to |x|^2 (Legendre polynomial of order 1, squared)
        print("# Using z-dipole distribution, random directions.")
      elseif _number_of_directions == 3 then
        -- Special case: simulate multlipole transition where angular density is proportional to square of (Legendre polynomial of order 3 or 4)
        print("# Using z-multipole distribution, random directions.")
      else
        print("# Using isotropic distribution, random directions.")
      end
    elseif _isotropic == 2 then --and group_index == 1 then -- print only for first group   
      print("# Using VMI-style distribution, both p_Z signs within each ray.")
    end
  end

  local local_ion_number = 1 + (ion_number-1) % _G.FLY2_IONS_PER_GROUP -- restarting from 1 in each new group if .fly2 file is set up to specify group size
  local group_index = 1 + math.floor((ion_number-1) / _G.FLY2_IONS_PER_GROUP)
  copy_of_ion_number = ion_number
  if group_index > highest_group_index then
	highest_group_index = group_index -- to let terminate_run() know how many groups there were (each group has its own statistics-variables)
  end
  
  if ion_number > _G.ions_flown_all_groups then -- To find out how many ions are flown (this adds the counts for different masses and energies)
    _G.ions_flown_all_groups = ion_number
  end
  if local_ion_number > _G.ions_flown then
    _G.ions_flown = local_ion_number
  end
  
  if local_ion_number == 1 then -- Only for the first ion (once or once per iteration if rerunning for different energy and/or mass)
    -- Reset the variables that track ion bunch hits
    if _number_of_directions < 1 then
      printf("_number_of_directions = %d is not allowed. Did you modify the wrong variable?", _number_of_directions);
      exit(2);
    end
    
    math.randomseed(random_seed) -- To make pseudorandom distributions generated in Lua reproducable (also same for each group in multi-group mode)
    -- NOTE: this doesn't affect the generator for SIMION FLY2 settings (e.g. Gaussian), which will differ between groups but seemingly be reproducable when re-flying
    -- So a _source_point_pattern>0 (e.g. 8) that uses Lua random number to get Gaussian should be used if reproducable results desired (e.g. for numeric optimization of lens parameters)
    
    --if #det_count < group_index then -- initialize arrays for this group (
	-- (Didn't find a way to do end-of-group printing of statistics, so now everything 
	--  is held in memory until segment.terminate_run which prints all of them)
    det_count[group_index] = {}
    det_y_min[group_index] = {}
    det_y_max[group_index] = {}
    det_y_mean[group_index] = {}
    det_y_Q[group_index] = {}
    det_y2_mean[group_index] = {}
    det_y2_Q[group_index] = {}
    det_t_mean[group_index] = {}
    det_t_Q[group_index] = {}
    for i = 1, _number_of_directions do
      det_count[group_index][i] = 0
      det_y_min[group_index][i] =  99999
      det_y_max[group_index][i] = -99999
      det_y_mean[group_index][i] = 0
      det_y_Q[group_index][i] = 0
      det_y2_mean[group_index][i] = 0
      det_y2_Q[group_index][i] = 0
      det_t_mean[group_index][i] = 0
      det_t_Q[group_index][i] = 0
    end
    
    
    if write_DLT == 1 or write_DLT == 4 then
      -- Start the output file
      output_DLT_file = io.open(output_filename, "wb");
      output_DLT_file:setvbuf("full") -- can select buffer size: , 16384)
      local head = "" -- Create rest of the the header in string form (to check its length)
      -- Start of aquisition timestamp (Ruby's ISO8601 doesn't include seconds for the timezone)
      head = head .. "\0" .. u24_as_string(2) --version
      head = head .. os.date("%Y-%m-%dT%H:%M:%S.000+00:00:00") -- time zone and milliseconds was not supported in this Lua
      -- Hardware settings:
      head = head .. "\4\255\1" -- Number of channels; Trigger channel (255=file was generated, not from hardware); Max hits per channel
      head = head .. "\0\0\0\0\0\0\0\0\0\0\0\0" -- Group range start, end; Trigger deadtime;
      head = head .. "\63\153\153\153\153\153\153\154" --Hardware time unit (DBL)
      head = head .. "\0" -- Property list (none)
      head = head .. "<Written by Lua-program in SIMION simulation>" --Arbitirary padding
      output_DLT_file:write("\61\30Psljus") --identifier
      output_DLT_file:write("\0" .. u24_as_string(8 + 4 + head:len()))
      output_DLT_file:write(head)
      
      _G.DLT_comment = ""
      _G.DLT_hit_count = 0
      _G.DLT_trigger_count = 0
    elseif write_DLT ~= 0 then
      -- Open the output file for resumed writing (multiple "Fly"s into one file)
      output_DLT_file = io.open(output_filename, "ab");
      output_DLT_file:setvbuf("full") -- can select buffer size: , 16384)
    else -- will not write to DLT file
      output_DLT_file = nil
    end
    
  end


  -- This is run for each particle.
  
  -- Set the start position:
  -- NOTE: Each source point is used _number_of_directions times (for each direction) before
  -- the next point is used, so to make sure all directions are represented in all points
  -- the number of particles needs to be a multiple of the number of directions.
  local point = floor((local_ion_number-1) / _number_of_directions) -- index to determine point position
  -- Select source point pattern:
  if source_point_pattern < 0 then
    -- Don't set the source point, use the position given by SIMION's FLY or FLY2 configuration.
    
  elseif source_point_pattern == 0 then
    -- Source points are spread only along r (SIMION: y). This was the only option before 2010-10-14.
    
    -- The number of points depends on number of particles which is not known, so this sequence is followed:
    -- (point,y_value/_source_point_spacing): (0,0), (1,1), (2,-1), (3,2), (4,-2), ...
    -- If the number of particles is an ODD multiple of _number_of_directions, the same number of points
    -- at y>0 and y<0 will be used.
    ion_py_mm = iif(point % 2 == 1, 1, -1) * ceil(point/2) * _source_point_spacing + source_point_y
    ion_px_mm = source_point_z
    
  elseif source_point_pattern == 5 or source_point_pattern == 7 then
    -- Source points are spread mainly along r (SIMION: y)
    -- but placed alternatingly at +- 0.05 mm  in z (SIMION: x) 
    -- The first, central point is actually at z=0 to make it symmetric, assuming odd number of points.
    -- When 11 points: one at z=0, five at z=-0.05, five at z=+0.05 the standard deviation is (by coincidence?) also 0.05.
    
    -- The number of points depends on number of particles which is not known, so this sequence is followed:
    -- (point,y_value/_source_point_spacing): (0,0), (1,1), (2,-1), (3,2), (4,-2), ...
    -- If the number of particles is an ODD multiple of _number_of_directions, the same number of points
    -- at y>0 and y<0 will be used.
    ion_py_mm = iif(point % 2 == 1, 1, -1) * ceil(point/2) * _source_point_spacing + source_point_y
    if point == 0 then -- the first
      ion_px_mm = source_point_z  -- center along z (SIMION: x)
    else
      if source_point_pattern == 7 then -- halfwidth 0.02 mm
        ion_px_mm = source_point_z + iif(floor(point/2)%2 == 0, 1, -1) * 0.02 --[mm]
      else -- halfwidth 0.05 mm: FWHM 0.1 mm
        ion_px_mm = source_point_z + iif(floor(point/2)%2 == 0, 1, -1) * 0.05 --[mm] % ORIGINAL
      end
    end
    
  elseif source_point_pattern == 1 then
    -- Source points are spread along both r (SIMION: y) and z (SIMION: x)
    
    -- Currently the number of particles is assumed to be 11*_number_of_directions, i.e. 11 points will be used.
    --y>0     A
    --       8 9
    --      4 3 5 
    --y=0: 1  0  2  ---> z
    --y<0    6 7  
    
    -- The spacing along z|x is 0.25 times the spacing along r|y
    local x_factor = 0.25
    if point < 6 then -- points 0 to 5
      ion_py_mm = source_point_y + floor(point / 3) * _source_point_spacing
      if point % 3 == 0 then -- point 0 and 3 are at source_point_z
        ion_px_mm = source_point_z
      else
        -- < or > source_point_z depending on point%3. Distance is 3 when floor(point/3)=0, distance is 2 otherwise.
        ion_px_mm = source_point_z + (-1)^(point%3) * (3 - floor(point/3)) * _source_point_spacing * x_factor
      end
    elseif point < 8 then -- points 6 & 7
      ion_py_mm = source_point_y - _source_point_spacing
      ion_px_mm = source_point_z + ((point%2)*2 - 1) * _source_point_spacing * x_factor
    elseif point < 10 then -- points 8 & 9
      ion_py_mm = source_point_y + 2 * _source_point_spacing
      ion_px_mm = source_point_z + ((point%2)*2 - 1) * _source_point_spacing * x_factor
    else -- point index 10 (the 11th point) and onwards: increase y, keep z (x) constant
      ion_py_mm = source_point_y + (3 + (point-10)) * _source_point_spacing
      ion_px_mm = source_point_z
    end
    
  elseif source_point_pattern == 2 then
    -- Source points are spread along both r (SIMION: y) and z (SIMION: x)
    
    --y<0:   C D
    --y>0:    B
    --y<0:   9 A
    --y>0:    8 
    --y<0:   6 7
    --y>0:  3 4 5 
    --y=0: 0  1  2  ---> z
    
    -- until 2013-10: The spacing along z|x is 0.25 times the spacing along r|y
    --local x_factor = 0.25
    -- until 2013-11-04: The spacing along z|x is 0.05 times the spacing along r|y. _source_point_spacing=1 approximately gives FWHM(z|x) 0.2 mm and range(y)=3mm
    local x_factor = 0.05
    if point < 6 then -- points 0 to 5
      ion_py_mm = source_point_y + floor(point / 3) * _source_point_spacing
      -- <, = or > source_point_z depending on point%3. Distance is 3 when floor(point/3)=0, distance is 2 otherwise.
      ion_px_mm = source_point_z + (point%3 - 1) * (3 - floor(point/3)) * _source_point_spacing * x_factor
    else
      if (point - 5) % 3 == 0 then -- point 8, 11, 14, ... : z-centered
        ion_px_mm = source_point_z
        ion_py_mm = source_point_y + ((point-5)/1.5 + 1) * _source_point_spacing
      else -- point 6, 7,  9, 10,  12, 13, ... : alternating low/high z-coordinate
        ion_px_mm = source_point_z + ((point%3)*2 - 1) * _source_point_spacing * x_factor
        ion_py_mm = source_point_y - (floor((point-5)/3)*2 + 2) * _source_point_spacing
      end
    end
  elseif source_point_pattern == 3 or source_point_pattern == 4 or source_point_pattern == 6 or source_point_pattern == 8 or source_point_pattern == 9 then
    -- Gaussian with std. _source_point_spacing along r (SIMION: y)
    -- 3: and std. 0.15 mm along z (SIMION: x)
    -- 4: and std. 0.02 mm along z (SIMION: x)), FWHM: 0.02 * (2*sqrt(2*log(2))) =approx= 0.0471 mm
    -- 8: and std. 0.084932 mm = 0.2mm / (2*sqrt(2*log(2))) along z (SIMION: x)
    -- 9: and std. 0.042466 mm = 0.1mm / (2*sqrt(2*log(2))) along z (SIMION: x), i.e. FWHM 100 um

    
    -- Using the algorithm in http://en.wikipedia.org/wiki/Box_Muller_transform
    local theta = TWO_PI * math.random()
    local R = math.sqrt(-2*math.log(1-math.random()))
    ion_py_mm = source_point_y + R * math.sin(theta) * _source_point_spacing
    if source_point_pattern == 4 then
      -- The standard deviation along z|x is fixed at 0.02 mm (FWHM ca 0.05 mm, i.e. 0.02*2.3548)
      ion_px_mm = source_point_z + R * math.cos(theta) * 0.02
    elseif source_point_pattern == 8 then
      -- The standard deviation along z|x is fixed at 0.084932 mm (FWHM 0.2 mm)
      ion_px_mm = source_point_z + R * math.cos(theta) * 0.084932
    elseif source_point_pattern == 9 then
      -- The standard deviation along z|x is fixed at 0.042466 mm (FWHM 0.1 mm)
      ion_px_mm = source_point_z + R * math.cos(theta) * 0.04246609 -- ORIGINAL
      --ion_px_mm = source_point_z + R * math.cos(theta) * 0.05/2.354820 -- given FWHM
    elseif source_point_pattern == 6 then
      -- The standard deviation along z|x is fixed at 0 mm (to compare how far 0.02 is from ideal)
      ion_px_mm = source_point_z
    else
      -- The standard deviation along z|x is fixed at 0.15 mm (seems very large now)
      ion_px_mm = source_point_z + R * math.cos(theta) * 0.15
    end
  else
    error("Unknown source pattern")
  end
  
  -- Speed
  if _G.CMD_kinetic_energy then
    -- Override GUI's or FLY2-file's energy with another value (from commandline.lua, or maybe only old version?)
    ion_mass = _G.CMD_mass
    ion_charge = _G.CMD_charge
    speed = ke_to_speed(_G.CMD_kinetic_energy, ion_mass)
  else
    -- Just get the speed that the particle has now, e.g. set by FLY2-file or GUI
    local tmp1, tmp2 
    speed, tmp1, tmp2 = rect3d_to_polar3d(ion_vx_mm, ion_vy_mm, ion_vz_mm)
  end
  _G.ion_energy[group_index] = speed_to_ke(speed, ion_mass) 
  if group_index <= #group_charge and group_charge[group_index] ~= ion_charge then
    error("Programming error makes it look like charge is not constant within particle group")
  end
  if group_index <= #group_mass and group_mass[group_index] ~= ion_mass then
    error("Programming error makes it look like mass is not constant within particle group")
  end
  group_charge[group_index] = ion_charge
  group_mass[group_index] = ion_mass
  group_charge[group_index] = ion_charge
  
  -- IMPROVEMENT: avoid hangup if source has been put outside simulation volume
  if (ion_py_mm-Y_CENTRE)^2+(ion_pz_mm-VERTICAL_CENTRE)^2 > R_GEOMETRY_MAX^2 or ion_px_mm <= 2 or ion_px_mm >= Z_GEOMETRY_MAX then
	print(string.format("FORCING SPLAT due to start outside: x(TOF) %.1f, y(optical) %.1f, z(vertical) %.1f mm. (optical centre %.1f, vertical centre %.1f mm)", ion_px_mm, ion_py_mm, ion_pz_mm, Y_CENTRE, VERTICAL_CENTRE))
    ion_splat = 1; --ion_py_mm = R_GEOMETRY_MAX*0.97; ion_px_mm = source_point_z;
  end
  
  -- Set the direction and vector velocity (using speed that was selected above)
  if _isotropic == 1 then
    -- Use isotropic (solid angle) direction distribution, not grouped by an enumeration of directions.
    -- Instead of using a deterministic count=sin(elevation)*max_count, a a a 3D-vector is picked
    -- inside a sphere, and normalized to the desired speed.
    
    -- TODO: can the resolution analysis be made by artificially binning the result into
    -- the usual set of (8 to 32) directions? Or some better method. Currently you need
    -- to write to DLT file and analyze it with the LabVIEW-programs as for experimental data.
    
    local a = math.random()-0.5; local b = math.random()-0.5; local c = math.random()-0.5; -- within a cube of half-side 0.5
    local r;
    if _number_of_directions == 2 then
	  -- Special case: simulate dipole or multipole transition where angular density is proportional
	  -- to the square of |Legendre polynomial| of some order.
	  
  	  if multipole_order == 0 then
		-- Special case to retain the (strange) implementation used 2016 (then there was no  multipole_order parameter):
		-- simulate dipole transition where angular density is proportional to |y|^2
		-- i.e. the spherical harmonic Y_1,0 (around y-axis insted of z).
		--  -- NOTE 2019: This version from 2016 has dipole along the laser propagation, not along the laser polarization which is z in SIMION (x is TOF).
		--  --            Was the choice of y made to keep y the main energy axis so that the long 
		--  --            Rayleigh range of light would be maximally challening for energy resolution?
		--  --            In a rotation-symmetric VMI it of course doesn't make much difference, but in the lab the polarization (and thus dipole) is perpendicular to the optical axis.
		-- Using the random number d=2*b in the range -1 to +1 to get the y-component of normalized direction vector,
		-- the probability density is y^2 <==> distribution function y^3
		-- <==> use the inverse function to pick y = sign(d) * abs(d)^(1/3) = sign(b) * abs(2*b)^(1/3) 
		  if b >= 0 then
			b =  (2*math.abs(b))^(1/3); -- directly get the y-component of 1-normalized direction vector
		  else
			b = -(2*math.abs(b))^(1/3); -- directly get the y-component of 1-normalized direction vector
		  end
		  -- Now only one more random number is needed to pick an angle in the (x,z)-plane (at the given y coordinate).
		  -- With b given, normalization gives 1^2 = b^2 + r_xz^2 where r_xz is radius in the (x,z)-plane
		  -- ==> r_xz = math.sqrt(1 - b^2).
		  -- Then x = r_xz * cos(uniformly_random_angle), z = r_xz * sin(same_uniformly_random_angle) ==>
		  -- Calculate the velocity
		  r = math.sqrt(1 - b^2);
		  ion_vx_mm = speed * r * math.cos(TWO_PI*a) -- angle from -pi to +pi;
		  ion_vy_mm = speed * b; -- along optical axis
		  ion_vz_mm = speed * r * math.sin(TWO_PI*a) -- angle from -pi to +pi;
		  
		  local tmp1, az, el
		  tmp1, az, el = rect3d_to_polar3d(ion_vy_mm, ion_vx_mm, ion_vz_mm)

	  elseif multipole_order >= 1 then
		-- NOTE 2019: The .lua program is now changed compared to the version used for making a poster-VMI-image and possibly optimizations in 2016.
		-- Then the y-axis was chosen, perhaps to have the energy spectrum with nonzero density along the optical axis (large source length) to make the focusing as difficult as possible?
		-- If such an orientation is desired, use  multipole_order = 0 instead of 1 now. See above.
		
	    -- Using the random number h_signed = 2*b = sign(b)*h in the range -1 to +1 to get the z-component of normalized direction vector,
		-- use the inverse function to pick b as z-componetn of 1-normalized direction vector, using approximation of C_3(h) or C_4(h).
		local h = 2*math.abs(b)
		  
		if multipole_order == 1 then
			-- Simulate dipole transition where angular density is proportional to |z|^2
			-- i.e. the spherical harmonic Y_1,0 (around z-axis).
			-- For a wavefunction with the Legendre polynomial of order 1 (Spherical harmonic with Y l=1, m=0) = z,
			-- the intensity is its modulus squared: |P_1(cos theta)|^2 = z^2.
			-- Using the random number d=2*b in the range -1 to +1 to get the z-component of normalized direction vector,
			-- the probability density is z^2 <==> distribution function y^3
			-- <==> use the inverse function to pick z = sign(d) * abs(d)^(1/3) = sign(b) * abs(2*b)^(1/3) 
			h =  (h)^(1/3); -- directly get the z-component of 1-normalized direction vector
			
		elseif multipole_order == 3 then
		  -- 2019-06: simulate multlipole transition where angular density 
		  -- is proportional to square of a Legendre polynomial of order 3 or 4 (not decided yet, will decide based on complex the expression and how the images look).
		  -- https://en.wikipedia.org/wiki/Legendre_polynomials :
		  --   I_3(z) = |P_3(cos theta)|^2 = |(5z^3 - 3z)/2|^2.
		  --   I_3(z) = |P_4(cos theta)|^2 = |(35z^4 - 30z^2 + 3)/8|^2.
		  -- Using Symbolic toolbox in Matlab: 
		  --   z = sym('z'); I_3 = ((5*z^3 - 3*z)/2)^2; I_4 = ((35*z^4 - 30*z^2 + 3)/8)^2; plot(subs(I_4, 'z',-1:0.01:1)); expand(I_4)
		  --                 I_7 = ((429*z^7 - 693*z^5 + 315*z^3 - 35*z)/16)^2;
		  -- The cumulative probability distributions are 
		  --   C_unnorm_3(z) = int(I_3,-1,'z');
		  --   C_unnorm_4(z) = int(I_4,-1,'z');
		  --   C_3(z) = C_unnorm_3(z) / C_unnorm_3(1) = (25*z^7)/8 - (21*z^5)/4 + (21*z^3)/8 + 1/2
		  --   C_4(z) = C_unnorm_4(z) / C_unnorm_4(1) = (1225*z^9)/128 - (675*z^7)/32 + (999*z^5)/64 - (135*z^3)/32 + (81*z)/128 + 1/2
		  -- Since there is no general inverse function for a polynomial, a search for articles suggests all
		  -- practical applications like this use a stored look-up table of the C_n(z) function (e.g. 1E5 points) and then interpolate linearly to get the inverse.
		  -- I don't need high accuracy and don't want to store a large array in Lua. 
		  -- a=sort([-1:0.02:1 -1.005:0.01:-0.85 1.005:-0.01:0.85 0.992 -0.992 0.982 -0.982 0.972 -0.972]); % points more dense near -1 and +1 to not give too much weight near 0.
		  -- C=double(subs(C_3,'z',a)); P_3 = polyfit(C-0.5,a,25);s=0:0.01:1; % order 23 is the first above 13 that is notably better, 25 further improved. Only odd orders are present/needed.
		  -- Podd=P_3; Podd(2:2:end)=0; P=Podd; plot(C,a,'.b', s,polyval(P,s-0.5),'-r'); axis([0 1 -1 1])
		  -- 	   P_3(1:2:end) = [2.71E+12, -4.96E+12, 4.01E+12, -1.89E+12, 5.73E+11, -1.17E+11, 1.63E+10, -1.55E+9 9.75E+7, -3.88E+6, 9.07E+4, -1.13E+3, 8.72]
		  --     are coefficients for z^25, z^23, ..., z^3, 1. 
		  --	Not final: h=0:0.001:0.5; plot(C,a,'.b', 0.5+h, (h).^(1/2.4) + 0.25*(h>0.28).*(h-0.28).^(1/5) + 0.06*exp(-((h-0.28)/0.2).^2).*(pi/2.2+atan((h-0.28)/0.02)),'-r'); axis([0.5 1 0 1])	
		  --	Not final: h=0:0.001:0.5; plot(C,a,'.b', 0.5+h, 0.8.*(0.8*h + 1.5*h.^2 + 4*h.^3 + 70*h.^5).^(1/3),'-r', 0.5+h, 0.750445 + 0.368*sign(h-0.278).*(h-0.278).^(1/3.875)); axis([0.5 1 0 1])
		  -- Nice approximation of right half of P_3 from the interval 0 to 1:
		  --	h = 0:0.001:1; plot(2*C-1,a,'.b', h, (h/3 + 0.15*h.^3 + 0.85*h.^4 - 0.4475*h.^7 + 0.133*(atan((h-0.555)/0.008)+pi/2).*exp(-(h-0.625).^2/0.33.^2) ).^(1/exp(1))); axis([0 1 0 1])
		  -- Decent approximation of right half of P_4 from the interval 0 to 1:
		  --  h=0:0.001:1; plot(2*C-1,a,'.b', h, 0.775*h - 0.15*h.^3 + 0.114*(atan((h-0.218)/0.009)+pi/2).*exp(-h/3.2) + 0.13*(atan((h-0.653)/0.0063)+pi/2).*exp(-h/0.787) ); axis([0 1 0 1])
	 
		  -- Using the random number h_signed = 2*b = sign(b)*h in the range -1 to +1 to get the z-component of normalized direction vector,
		  -- use the inverse function to pick b as z-componetn of 1-normalized direction vector, using approximation of C_3(h) or C_4(h).
		  
		  -- For third order Legendre polynomial: Generate random number b with distribution like |P_3|^2 
		  h = (h/3 + 0.15*h^3 + 0.85*h^4 - 0.4475*h^7 + 0.133*(math.atan((h-0.555)/0.008)+TWO_PI/4)*math.exp(-(h-0.625)^2/0.33^2) )^(1/math.exp(1));
		  if h > 1 then 
		    -- The approximative formula can give h slightly larger than 1 (gives error in math.sqrt below). Truncate it to 1.
			h = 1;
		  end
		elseif multipole_order == 4 then
		  -- For fourth order Legendre polynomial: Generate random number b with distribution like |P_4|^2 
		  h = 0.775*h - 0.15*h^3 + 0.114*(math.atan((h-0.218)/0.009)+TWO_PI/4)*math.exp(-h/3.2) + 0.13*(math.atan((h-0.653)/0.0063)+TWO_PI/4)*math.exp(-h/0.787);
		  if h > 1 then 
		    -- The approximative formula can give h slightly larger than 1 (gives error in math.sqrt below). Truncate it to 1.
			h = 1;
		  end
		elseif multipole_order == 7 then
		  -- Giving up analytic approximations, using lookup table with 17 points for C_7:
		  --   h = 0:0.0625:1; interp_7 = interp1(2*C-1,a,h); plot(2*C-1,a,'-b', h, interp_7,'.r-' ); axis([0 1 0 1]),h(end)
		  local interp = {0, 0.149, 0.2015, 0.2513, 0.3235, 0.524, 0.5715, 0.6121, 0.6609, 0.8132, 0.8486, 0.8743, 0.901, 0.97517, 0.98809, 0.99521, 1}
		  local w, high
		  high = math.ceil(h / 0.0625)
		  if high <= 0 then
		    h = 0;
		  elseif high > 16 then
		    h = 1;
		  else -- linear interpolation
		    w = high - h / 0.0625
			h = w*interp[high] + (1-w)*interp[1+high] -- 1-based array indexing in Lua
		  end
		  
        else
		  printf("multipole_order %d is not supported", multipole_order)
		  error("The given value for multipole_order is not supported.")
	    end
		
		if b >= 0 then
		  b = h
		else
		  b = -h
		end
		-- Now only one more random number is needed to pick an angle in the (x,y)-plane (at the given z coordinate).
		-- With b given, normalization gives 1^2 = b^2 + r_xy^2 where r_xy is radius in the (x,y)-plane
		-- ==> r_xy = math.sqrt(1 - b^2).
		-- Then x = r_xy * cos(uniformly_random_angle), z = r_xy * sin(same_uniformly_random_angle) ==>
		-- Calculate the velocity
		r = math.sqrt(1 - b^2);
		ion_vx_mm = speed * r * math.cos(TWO_PI*a) -- angle from -pi to +pi;
		ion_vy_mm = speed * r * math.sin(TWO_PI*a) -- angle from -pi to +pi;
		ion_vz_mm = speed * b; -- along the axis that is neither optical nor TOF, i.e. the vertical axis in the lab
	  end
		
    else -- Normal case: fully isotropic
      r = a*a + b*b + c*c
      while r < 1E-6 or r > 0.25 do
        -- When not inside the ball of radius 0.5 (or when r is very small): pick another random vector
        a = math.random()-0.5; b = math.random()-0.5; c = math.random()-0.5; -- within a cube of half-side 0.5
        r = a*a + b*b + c*c
      end
      -- Now [a,b,c] is an unnormalized randomized direction vector
      -- Calculate the velocity
      r = speed / math.sqrt(r) --normalization factor (the sqrt was delayed until here, to not repeat it if the loop runs several iterations)
      ion_vx_mm = a * r; ion_vy_mm = b * r; ion_vz_mm = c * r;
    end

    if _number_of_directions > 2 then -- log each direction only when grouping into different directions will be done at end
      -- Round the angle to the bin index nearest this direction
      local tmp1, az, el
      tmp1, az, el = rect3d_to_polar3d(ion_vy_mm, ion_vx_mm, ion_vz_mm)
      -- The second coordinate is along the axis described by the returned el
      -- The above doesn't give coordinates the way I want them, it uses el = +-90 on the poles and 0 on the equator.
      --   (vx=1,vy=0,vz=0) (to detector)    ==> el = 90, az =   0,  want 0
      --   (vx=-1,vy=0,vz=0) (away fron det.)==> el =-90, az =   0,  want 180
      --   (vx=0,vy=1,vz=0) (just sideways)  ==> el =  0, az =   0,  want 90
      --   (vx=0,vy=-1,vz=0) (just sideways) ==> el =  0, az =-180,  want 270
      --   (vx=0,vy=0,vz=1)                  ==> el =  0, az = -90,  want 90 or 270
      --   (vx=2,vy=0,vz=1)                  ==> el = 63, az = -90,  want 27 or 333
      --   (vx=-1,vy=0.7,vz=0.7)             ==> el =-45, az =  45,  want 225
      -- To comply with the binning used when _isotropic==0, I want to ignore "az" and use a span from 0 to 360 for el
      if az >= 0 then
        el = 90 - el
      else
        el = 270 + el
      end
      if _G.ion_direction_index == nil then
        _G.ion_direction_index = {}
        _G.ion_direction_unit_y = {}
        _G.ion_direction_unit_z = {}
      end
      -- The center angle for a bin with index i is (360/_number_of_directions) * (i+1).
      _G.ion_direction_index[local_ion_number] = math.floor(0.5 + ((el/360)*_number_of_directions)) % _number_of_directions + 1
    
      tmp1 = sqrt(b*b + c*c) -- Save initial velocity's unit vector in plane parallell to detector
      _G.ion_direction_unit_y[local_ion_number] = b/tmp1
      _G.ion_direction_unit_z[local_ion_number] = c/tmp1
      --printf("e %.0f \t a %.0f \t i %d \t ion_vx_mm %.1f, ion_vy_mm %.1f", el, az, _G.ion_direction_index[local_ion_number], ion_vx_mm, ion_vy_mm)
    end
  
  else -- (supported values of _isotropic: 0, -1, 2, 3, 4)
    -- A discrete set of elevation-directions. This is the old & usual way.
    local azimuth = 0;

    if _isotropic == 0 then
      -- Elevations in full (R,Z)-momentum plane. For 3D momentum imaging (unlike 2D VMI).
      -- For _number_of_directions=8 it gives 0 45 90 135 180 225 270 315 [degrees].
      -- _number_of_directions is typically a multiple of four.
      elevation = (360/_number_of_directions) * ((local_ion_number-1) % _number_of_directions) --[degrees]

    elseif _isotropic == 2 then
      -- Elevations grouped without sign of Z-momentum. For 2D VMI (unlike 3D).
      -- Typically used with even _number_of_directions (e.g. 8) to avoid elevation=0 where VMI info is not possible.
      -- Within each direction, every second ion is directed with negative p_Z (away from detector) and every second with positive p_Z.
      -- For the directions +90 and -90, the alternating p_Z signs have no effect (p_Z=0).
      
      -- Distinct rays for p_Y>0 and p_Y<0 (to be able to test y-offset sensitivity easily):
      local rel_elev = (math.floor((local_ion_number-1)/2) % _number_of_directions) * (180/(_number_of_directions-1)) --[degrees]
      -- Alternative, only at p_Y>0 (to not waste simulation effort in symmetric case):
      --rel_elev = (math.floor((local_ion_number-1)/2) % _number_of_directions) * (90/(_number_of_directions-1)) --[degrees]

      if (local_ion_number % 2) == 1 then -- alternate between positive and negative p_Z sign
        elevation = 90 - rel_elev;
      else
        elevation = 90 + rel_elev;
      end

      -- NOTE: the 90 degree direction is its own mirror-image in the p_Z=0 plane, but two particles are run
      -- for that direction anyway. This ensures that standard deviation etc. of this ray (direction)
      -- is comparable to standard deviations for other directions (deviations represent only effects
      -- of lens and source point distribution). If particles are treated without grouping by rays,
      -- it may also be seen as an increased weight for the 90-degree direction which is reasonable
      -- when the polarization is perpendicular to the spectrometer axis (as is the case for 2D VMI).
    
    elseif _isotropic == 3 or _isotropic == 4 then
      -- Similar to the VMI _isotropic==2, but the number of particles per elevation-value is variable
      -- to approximate the sin(elevation) dependence to get isotropic solid-angle distribution.
      -- Elevation is 0 to 180 degrees (symmetric around 90), azimuth is spanning 0 to 360 degrees.
      -- Both angles are set deterministically as function of local_ion_number.
      -- To ensure VMI resolution can be calculated for all directions (i.e. exclude 0 degrees), _number_of_directions should be even.
      -- _isotropic == 3: First elevation is 90 degrees. ==> direction_elevation[1] == 180 - direction_elevation[1], i.e. p_Z=0 and this direction will not be mirrored.
      -- _isotropic == 4: First elevation is lower than 90 degrees. The average of direction_elevation[1] and its mirror 180-direction_elevation[1] is 90 degrees.
      -- _isotropic == 3 is recommended since it ensures the "worst case" (theta=+-90) for transverse velocity is tested and it is more similar to _isotropic == 2.
      
      if ion_number == 1 then -- First ion of all
        -- Prepare array with the desired count for each direction (since it varies between directions here)
        -- and also store the elevation (theta) angle of each direction (to not need re-computations).
        direction_elevation = {};
        local share_per_direction = {}; -- proportional to ideal number of ions for the direction, before rounding
        ions_per_direction = {}; -- after rounding with adjustments to reach desired total count
        local quot = {}; -- temporary quantity used in Sainte-Laguë's method

        if _number_of_directions < 2 then
          error("Can not use quasi-isotropic mode with discrete polar angles (_isotropic == 3 or 4) with _number_of_directions < 2.");
        end

        if _isotropic == 4 then -- variant
          -- +-90 degrees is not reached by the highest elevation, meaning all elevations can have alternating p_Z signs
          for dir = 1, _number_of_directions do
            direction_elevation[dir] = 90 - (dir - 0.5) * 180/_number_of_directions;
            -- Since each direction is used with alternating signs for p_Z (forward and backward), the count is doubled after rounding:
            share_per_direction[dir] = 2 * math.abs(math.cos((90 - (dir-1)*180/_number_of_directions)*DEGREE) - math.cos((90 - dir*180/_number_of_directions)*DEGREE)) / 2;
          end
        else -- main option
          -- +-90 degrees is the highest elevation used, making it the only elevation where p_Z is zero (not alternating sign).
          for dir = 1, _number_of_directions do
            direction_elevation[dir] = 90 - (dir - 1) * 180/(_number_of_directions-1);
            -- Most directions are used with alternating signs for p_Z (forward and backward), doublign the count after rounding:
            share_per_direction[dir] = 2 * math.abs(math.cos((90 - (dir-0.5)*180/(_number_of_directions-1))*DEGREE) - math.cos((90 - (dir-1.5)*180/(_number_of_directions-1))*DEGREE)) / 2;
            if dir == 1 or dir == _number_of_directions then
              -- The +90 and -90 degrees are not doubled by alternating p_Z sign, thus the above expression did not need the doubling.
              share_per_direction[dir] = share_per_direction[dir] / 2;
            end
          end
        end
        
        -- Now need to adjust the rounding of ions_per_direction[dir] to ensure sum is _G.FLY2_IONS_PER_GROUP, otherwise group index determination doesn't work!
        -- Use the "Sainte-Laguë method" ["uddatalsmetoden" in Swedish] to distribute the share onto an integer number of total ions.
        -- Initialization
        local remaining = _G.FLY2_IONS_PER_GROUP;
        for dir = 1, _number_of_directions do
          if direction_elevation[dir] == -90 or direction_elevation[dir] == 90 then
            ions_per_direction[dir] = 1; -- ensure each direction gets at least one particle
          else -- for directions where both p_Z>0 and p_Z<0 should be used
            ions_per_direction[dir] = 2; -- ensure each direction gets at least two particles
          end
          quot[dir] = share_per_direction[dir] / (2*ions_per_direction[dir] + 1);
          remaining = remaining - ions_per_direction[dir];
        end
        if remaining < 0 then
          error(string.format("With _number_of_directions=%d, %d particles is not sufficient to ensure proportionality.", _number_of_directions, _G.FLY2_IONS_PER_GROUP));
        end
        while remaining > 0 do
          -- For small arrays, math.max(unpack(the_array)) can be used to find maximum. Tested that max number of directions where it works is at least 17 in SIMION.
          local m = math.max(unpack(quot));
          local d = 0; -- Get the first (or last if remaining is even) index among those with the highest quotient.
          if (remaining % 2) == 1 then
            for i = 1, _number_of_directions do
              if quot[i] >= m then
                d = i; break;
              end
            end
          else
            for i = _number_of_directions, 1, -1 do -- reversed iteration
              if quot[i] >= m then
                d = i; break;
              end
            end
          end
          if d == 0 then
            error("Failed to find which index has the maximum quotient %.9g", m);
          end
          -- Now we have chosen to award one more particle to the direction with index d.
          ions_per_direction[d] = ions_per_direction[d] + 1;
          remaining = remaining - 1;
          -- Update its quotient for the next round.
          quot[d] = share_per_direction[d] / (2*ions_per_direction[d] + 1);
        end
        -- Now the desired number of particles have been distributed among the directions!
        
        local msg = ""
        for dir = 1, _number_of_directions do
          --printf("# Direction %d is %.1f degrees with (combining all p_Z signs) %d particles.", dir, direction_elevation[dir], ions_per_direction[dir]);
          if dir == 1 then
            msg = string.format("%d", ions_per_direction[dir]);
          else
            msg = string.format("%s,%d", msg, ions_per_direction[dir]);
          end
        end
        printf("# Using VMI-style 3D-isotropic, both p_Z signs within each ray. Count per ray: [%s]", msg);
      end
      if _G.ion_direction_index == nil then
        _G.ion_direction_index = {};
        _G.ion_direction_unit_y = {};
        _G.ion_direction_unit_z = {};
      end

      -- Choose direction by stepping through the directions and comparing their desired count with the (local) ion index.
      local dir_index = nil;
      local last_ion_for_current_direction = 0;
      for dir = 1, _number_of_directions do
        last_ion_for_current_direction = last_ion_for_current_direction + ions_per_direction[dir];
        if local_ion_number <= last_ion_for_current_direction then
          dir_index = dir;
          break;
        end
      end
      _G.ion_direction_index[local_ion_number] = dir_index;
      local number_within_direction = local_ion_number - (last_ion_for_current_direction-ions_per_direction[dir_index]); -- from 1 to ions_per_direction[dir]
      -- Alternation between positive and negative p_Z sign (has no effect if direction_elevation[dir] == 90)
      if (local_ion_number % 2) == 1 then
        elevation = direction_elevation[dir_index] % 360;       -- ranges from 0 to +90 and 270 to 360 degrees (p_Z >= 0)
         azimuth = (((number_within_direction-1 + local_ion_number/33) / max(3, ions_per_direction[dir_index])) % 1) * 360 ;
      else
        elevation = 180 - direction_elevation[dir_index]; -- ranges from  90 to 270 degrees (p_Z <= 0)
        azimuth = (((number_within_direction-1 + local_ion_number/73) / max(3, ions_per_direction[dir_index]) + 0.375) % 1) * 360 ; -- 0.625=225/360 offset to use other angles than in forward direction
      end
      -- Alternative patterns for the azimuth (discretized uniform distribution)
      -- azimuth = ((local_ion_number / max(3, ions_per_direction[dir_index])) % 1) * 360 ;
      -- azimuth = 360 * math.random(); --[degrees] randomly, uniformly distributed -- random

      --printf("# ion #%d direction index %d : theta %.2f degrees, phi %.2f degrees", local_ion_number, dir_index, elevation, azimuth); --DEBUG
        
      if elevation < 180 then
        _G.ion_direction_unit_y[local_ion_number] = math.cos(azimuth*DEGREE)
        _G.ion_direction_unit_z[local_ion_number] = -math.sin(azimuth*DEGREE)
      else
        _G.ion_direction_unit_y[local_ion_number] = -math.cos(azimuth*DEGREE)
        _G.ion_direction_unit_z[local_ion_number] = math.sin(azimuth*DEGREE)
      end
      
    elseif _isotropic == -1 then
      -- Just randomize the angle around the spectrometer axis. This gives too
      -- high solid angle density at the "poles" (elevation 0 and 180 degrees)
      -- compared with an isotropic solid-angle distribution. The same over-density at poles
      -- is actually present also when _isotropic==0 so that azimuth is kept constnt.
      -- Note the range for azimuth is 180 degrees, since the elevation has a 360 degree range here.
      azimuth = 180 * math.random() --[degrees]

      -- Save initial velocity's unit vector in plane parallell to detector
      if _G.ion_direction_unit_y == nil then
        _G.ion_direction_unit_y = {}
        _G.ion_direction_unit_z = {}
      end
      if elevation < 180 then
        _G.ion_direction_unit_y[local_ion_number] = math.cos(azimuth*DEGREE)
        _G.ion_direction_unit_z[local_ion_number] = -math.sin(azimuth*DEGREE)
      else
        _G.ion_direction_unit_y[local_ion_number] = -math.cos(azimuth*DEGREE)
        _G.ion_direction_unit_z[local_ion_number] = math.sin(azimuth*DEGREE)
      end
    end

    -- Old: ion_vx_mm, ion_vy_mm, ion_vz_mm = polar3d_to_rect3d(speed, 0, elevation)
    -- Calculate the velocity
    -- NOTE order: the elevation is the "theta" and azimuth is the "phi" (at theta=0 or 180, the azimuth doesn't matter)
    ion_vy_mm, ion_vx_mm, ion_vz_mm = polar3d_to_rect3d(speed, azimuth, 90 - elevation)
    --printf("#%d : V (%.1f,%.1f,%.1f) : D (,%.1f,%.1f)", local_ion_number, ion_vx_mm, ion_vy_mm, ion_vz_mm, _G.ion_direction_unit_y[local_ion_number], _G.ion_direction_unit_z[local_ion_number]); --DEBUG
    
  end
  
  -- Initial offset velocity, to simulate skimmed pulsed supersonic jet (NOTE: vertical molecular beam direction is z in current SIMON geometries, could be called "x" in the swapped coordinates used for many variables)
  -- 2019-10-08 Apparently (in the z-unsymmetric both8ue.iob) a POSITIVE ion_vz_mm makes the velocity go down (towards lower Z). Is there any Lua-rotation affecting it?
  ion_vy_mm = ion_vy_mm + initial_vy -- NOTE the unit [mm/us] = [km/s] 
  ion_vz_mm = ion_vz_mm + initial_v_vertical -- NOTE the unit [mm/us] = [km/s]
  
  --print(string.format("ion %d, pg %d, tot %d, group %d", local_ion_number, _G.ions_flown, _G.ions_flown_all_groups, group_index)) --DEBUG
end

function segment.init_p_values()
  if #det_count == 0 then
    error("ERROR: segment.init_p_values() was called before segment.initialize(). This happens if the particles start outside or at edge of workbench, i.e. error in .fly2-file or source_point_z.")
  end

  if DoNotSetPotentials == 0 or DoNotSetPotentials == -1 then
    -- Assign the values to the Simion potential array (which we have already defined in the U[])
    -- When DoNotSetPotentials == -1 the _U_... variables are not used, just assuming
    -- external script (adjust.lua) has set the U-variable
    local start_index = 1
    if WALLS > 1 then
      -- There is no electrode at index 1 for walls
      start_index = WALLS;
    end
    for i = start_index, #U do
      if U[i] ~= nil and adj_elect[i] ~= nil then
        adj_elect[i] = U[i]
      end
    end
    potential_display_needs_updating = 1 -- request update of potential display
    
  else
    if ion_instance == 1 then
      -- Print message only once (not repeated if multiple potential array instances in workbench)
      print("# Not setting electrode potentials, using fast-adjusted values.")
      printf("U_AMCP = %g, U_Alast = %g, UBfirst = %g", adj_elect[A_MCP] or 0, adj_elect[A_LAST] or 0, adj_elect[B_FIRST] or 0)
    end
    local electrode_count = math.max(A_LAST, math.max(A_MCP, B_MCP)); --NOTE: needs to know the number of electrodes
    for i = 1, electrode_count do -- More informative version:
      -- NOTE this will be run once per potential array instance, but each instance doesn't have all the electrode indices
      -- and there may be some indexes not used at all.
      if adj_elect[i] ~= nil then
        
        U[i] = adj_elect[i]; -- read the user-adjusted potentials instead, to define U[] (used by mass-resolution calculator)
        if i < electrode_count and adj_elect[i+1] ~= nil then
          --OLD: printf("U[%d] = %g,  V[%d] = %g", i, adj_elect[i+1], i, adj_elect[i+1] - adj_elect[i+1+1])
          printf("U[%d] = %g,  V[%d] = %g", i, adj_elect[i], i, adj_elect[i] - adj_elect[i+1])
        else
          --OLD: printf("U[%d] = %g,  (is last)", i, adj_elect[i+1], i)
          printf("U[%d] = %g,  (no next)", i, adj_elect[i], i)
        end
      end
    end
    
  end
end

-- Update PE surface display.
function segment.other_actions()
  
  if potential_display_needs_updating ~= 0 then
    potential_display_needs_updating = 0
    sim_update_pe_surface = 1           -- update the PE surface display (sadly not always the curve view)
    redraw_screen() -- redraw the current view (e.g. potential curves)
  end

  -- Check if ion has splatted
  if ion_splat ~= 0 then
    --(for cylinder symm B-side) if abs(ion_px_mm - detectorA_x) < 1 or (ion_px_mm - detectorB_x) < 16 then -- hit the detector
    --For more general B-side (2D slope allowed)
    if abs(ion_px_mm - detectorA_x) < 1 or (abs(detectorB_linecoeff[1]*ion_px_mm + detectorB_linecoeff[2]*ion_py_mm + detectorB_linecoeff[3]) <= 2.5) then -- hit the detector
      local det_centre_y = 0;
      if abs(ion_px_mm - detectorA_x) >= 1 and #detectorB_linecoeff >= 4 then
        -- Since not near detector A, the particle must have hit detector B (which may have nonzero y-centre)
        det_centre_y = detectorB_linecoeff[4];
      else
        det_centre_y = detectorA_y;
      end
      
      local local_ion_number = 1 + (ion_number-1) % _G.FLY2_IONS_PER_GROUP -- restarting from 1 in each new group if .fly2 file is set up to specify group size
      local group_index = 1 + math.floor((ion_number-1) / _G.FLY2_IONS_PER_GROUP)
      
      if _isotropic == 1 or _isotropic == 3 or _isotropic == 4 then -- random or deterministic, (approximately) isotropic directions. The initial direction has been classified into to the most appropriate bin.
        if _isotropic == 1 and _number_of_directions <= 2 then
          -- For _number_of_directions==1 and also for _number_of_directions==2 (which are used for Legendre-polynomial-weighted isotropy (dipole or multipole transition)),
          -- the direction is not interesting (and speed up when simulating VMI image with many particles).
          dir = 1; -- don't care about direction, use ray 1 always
          -- NOTE: for users who want to log impact positions (projected onto _G.ion_direction_unit...) while _isotropic=1, use _number_of_directions>2 (e.g. 4).
        else
          if _G.ion_direction_index == nil or _G.ion_direction_index[local_ion_number] == nil then
            error(string.format("ERROR: missing _G.ion_direction_index[%d] whith _isotropic=1.", local_ion_number))
          else
            dir = _G.ion_direction_index[local_ion_number];
          end
        end
      elseif _isotropic == 2 then
        -- Deterministic elevations grouped without sign of Z-momentum. For 2D VMI (unlike 3D).
        dir = math.floor((local_ion_number-1)/2) % _number_of_directions + 1 -- direction index (starting at 1 here)
      else -- (0 or -1) Deterministic elevations from full (Z,R)-momentum plane.
        dir = ((local_ion_number-1) % _number_of_directions) + 1 -- direction index (starting at 1 here)
      end

      local position -- the radial or y coordinate used for VMI resolution assessment
      if _isotropic == 1 or _isotropic == -1 or _isotropic == 3 or _isotropic == 4 then
        if _number_of_directions > 2 then -- only if we care about directions
          -- Particles may have direction (initial velocity) along not just along y but also SIMION-z-axis (the other axis perpendicular to flight axis). 
          -- Project impact position to an axis parallell to initial velocity, and in detector plane, to get the "radius" relevant for VMI.
          -- For an initial velocity with (v0_y>0, v0_z=0) this is the usual y-axis projection.
          -- An "overfocused" trajectory may still end up with a negative coordinate value, meaning "on wrong side of detector middle for VMI".
          position = _G.ion_direction_unit_y[local_ion_number] * (ion_py_mm - det_centre_y)
                   + _G.ion_direction_unit_z[local_ion_number] * (ion_pz_mm - VERTICAL_CENTRE)
          -- TODO: this still doesn't account for possible slope of detector, which is encoded in detectorB_linecoeff...
          
          --printf("#%d dir%d: (%.1f,%.1f,%.1f) --> R %.2f", local_ion_number, dir, ion_px_mm, ion_py_mm, ion_pz_mm, position); --DEBUG
          --printf("Relative %.3f", position / math.sqrt((ion_py_mm-det_centre_y)^2 + (ion_pz_mm-VERTICAL_CENTRE)^2)) --is 1.000 --DEBUG
        else
          position = 0; -- ignoring position information in Lua, may analyse full log file from Matlab instead
        end
      else --: Particles have initial velocity only along +- y axis. Use signed y-coordinate as position
        position = ion_py_mm - det_centre_y
        -- NOTE: to be fully equivalent with the _isotropic ~= 0 case, should switch sign when
        -- dir > _number_of_directions/2, but I don't think any resolution-value can differ (only sign on y-mean) so for performance I don't do that.
      end
      -- NOTE: when _source_point_spacing ~= 0 (and the point pattern is not rotation-symmetric in y,z-plane)
      -- the results from _isotropic=-1 and isotropic = 0 differ slightly. Otherwise they are equal.
      -- The results with _isotropic=1 differ in many ways, probably just because the isotropic distribution
      -- gives different weight (particle count) for the different directions.

      det_count[group_index][dir] = 1 + det_count[group_index][dir]
      if position < det_y_min[group_index][dir] then
        det_y_min[group_index][dir] = position
      end
      if position > det_y_max[group_index][dir] then
        det_y_max[group_index][dir] = position
      end

-- Compute standard deviation in one pass, in a way that attempts to avoid numeric errors
-- (difference between large similar numbers). Reference: Donald E. Knuth (1998). "The Art
-- of Computer Programming, volume 2: Seminumerical Algorithms, 3rd edn.", p. 232. Addison-Wesley,
-- via http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance--On-line_algorithm
      -- Update mean and variance for detected y position
      delta = position - det_y_mean[group_index][dir]
      det_y_mean[group_index][dir] = det_y_mean[group_index][dir] + delta/det_count[group_index][dir]
      det_y_Q[group_index][dir] = det_y_Q[group_index][dir] + delta*(position - det_y_mean[group_index][dir]) --Yes, the new value of _mean needs to be used with the old delta
      -- Update mean and variance for square of detected y position (proportional to kinetic energy if p_Z=0)
      delta = position*position - det_y2_mean[group_index][dir]
      det_y2_mean[group_index][dir] = det_y2_mean[group_index][dir] + delta/det_count[group_index][dir]
      det_y2_Q[group_index][dir] = det_y2_Q[group_index][dir] + delta*(position*position - det_y2_mean[group_index][dir]) --Yes, the new value of _mean needs to be used with the old delta
      -- Update mean and variance for detected time of flight
      delta = ion_time_of_flight - det_t_mean[group_index][dir]
      det_t_mean[group_index][dir] = det_t_mean[group_index][dir] + delta/det_count[group_index][dir]
      det_t_Q[group_index][dir] = det_t_Q[group_index][dir] + delta*(ion_time_of_flight - det_t_mean[group_index][dir]) --Yes, the new value of _mean needs to be used with the old delta

      if write_DLT ~= 0 then
        -- Write the hit to DLT file, using format version 2.
        -- Since hits are written at once, there is currently no way to use hit index>0 (could have used one index per direction, per source point or such)
        -- Compute the "hardware times" for each channel
        -- TODO: check which sign etc. is used, so that the correct x and y are obtained in the program
        local t_x = (ion_pz_mm-VERTICAL_CENTRE) / G_x
        local t_y = (ion_py_mm-detectorA_y) / G_y
        local t_z = ion_time_of_flight + time_offset_DLT
        output_DLT_file:write("\255\0\0\4") -- the event(group) marker and size
        output_DLT_file:write("\127\255\255\255\255\255\255\255") -- NaN absolute trigger time
        output_DLT_file:write("\0" .. u24_as_string(math.floor(0.5
                            + (t_z - t_x)/HW_TICK) )) --time for channel 0
        output_DLT_file:write("\1" .. u24_as_string(math.floor(0.5
                            + (t_z + t_x)/HW_TICK) )) --time for channel 1
        output_DLT_file:write("\2" .. u24_as_string(math.floor(0.5
                            + (t_z - t_y)/HW_TICK) )) --time for channel 2
        output_DLT_file:write("\3" .. u24_as_string(math.floor(0.5
                            + (t_z + t_y)/HW_TICK) )) --time for channel 3
        _G.DLT_hit_count = _G.DLT_hit_count + 1 -- count actual number of written hits
      end
    
    else -- missed the detector
      --print("MISS by -- at r=-- d=--", local_ion_number, ion_py_mm, ion_px_mm)
    end
  else
    -- NOTE: if a particle exits the workbench geometry without splatting (no electrode surrounding all of it)
    -- the .terminate() segment will not be called. If this happens for last particle in group,
    -- the group summary won't be printed!
    --if ion_px_mm >= Z_GEOMETRY_MAX or ion_px_mm <= 1 or math.abs(ion_py_mm) > R_GEOMETRY_MAX then -- not checking math.abs(ion_pz_mm) > R_GEOMETRY_MAX, since this dimension is typically not used
	if ion_px_mm >= Z_GEOMETRY_MAX or ion_px_mm <= 1 or (ion_py_mm-Y_CENTRE)^2+(ion_pz_mm-VERTICAL_CENTRE)^2 > R_GEOMETRY_MAX^2 then
      -- Force a splat when particle is nearly outside simulation volume
      ion_splat = 1
	  --print(string.format("FORCING SPLAT for particle leaving volume: x(TOF) %.1f, y(optical) %.1f, z(vertical) %.1f mm. (optical centre %.1f, vertical centre %.1f mm)", ion_px_mm, ion_py_mm, ion_pz_mm, Y_CENTRE, VERTICAL_CENTRE))
    end
    
  end
end

function segment.terminate_run()
  -- When finishing the last run.
  -- It was not reliable enough to use segment.terminate() for the last particle to print,
  -- because in an open geometry the last particle may escape without hitting anything.
  -- Newer SIMION has terminate_run() which allows us to print all the group summaries 
  -- in a loop at the end (the statistics is kept separate by indexing by group_index).
  for g = 1, highest_group_index do 
    print_summary(g, g == highest_group_index)
  end
  
  -- To aid Matlab analysis (needing to know where centre of detectors are to analyze trajectories),
  -- print also the detector centeres (and maybe more)
  print(string.format("Detector coordinates: A %.1f, %.1f, %.1f; B %.2f, %.4f, %.4f, %.4f, %.2f;", detectorA_x, detectorA_y, VERTICAL_CENTRE, detectorB_x,
		detectorB_linecoeff[1] or 0, detectorB_linecoeff[2] or 0, detectorB_linecoeff[3] or 0, detectorB_linecoeff[4] or 0))
end

function segment.terminate()
-- NOTE: if a particle exits the workbench geometry without splatting (no electrode surrounding all of it)
-- the .terminate() segment will not be called. If this happens for last particle in group,
-- the group summary won't be printed!

  -- To prevent Simion from restoring the old potential values when the flying is completed
  sim_retain_changed_potentials = 1
  
  --local local_ion_number = 1 + (ion_number-1) % _G.FLY2_IONS_PER_GROUP -- restarting from 1 in each new group if .fly2 file is set up to specify group size
  --local group_index = 1 + math.floor((ion_number-1) / _G.FLY2_IONS_PER_GROUP)
  --print(string.format("DEBUG: Terminating g %d, l %d, i %d, copy %d", group_index, local_ion_number, ion_number, copy_of_ion_number)) --DEBUG
  
  --if ion_number == _G.ions_flown or local_ion_number == _G.FLY2_IONS_PER_GROUP then -- when terminating for the last ion (or last in group)
    -- (If grouping, then ions_flown is the same as _G.FLY2_IONS_PER_GROUP)
	--print_summary(group_index, ion_number == _G.ions_flown_all_groups) -- not needed anymore, now done from segment.terminate_run()
  --end
  
  --NOTE: the file writing could be put here instead, possibly with the problem (is that why "other_action" was used for calculations above?)
  --that maybe there is no way to close the file after the last particle has splatted.
  --print(string.format("Terminated at %.3f, %.3f, %.3f", ion_px_mm, ion_py_mm, 1000*ion_time_of_flight))
end

function print_summary(group_index, is_last)
  --print(string.format("DEBUG: print_summary g %d, i-copy %d. m %.3g, ch %.3g", group_index, copy_of_ion_number, group_mass[group_index], group_charge[group_index]))
    
    local points_used = _G.ions_flown / _number_of_directions; -- (average) number of particles per direction
    
    local flown_message
    if source_point_pattern < 0 and _G.FLY2_SOURCE_DISTR ~= nil then
      -- When using Gaussian source distribution from .fly2 file, a string gives its sizes on the form "[std_axial,std_SIMION_Y,std_SIMION_Z]"
      flown_message = string.format("Flew %d ions @ %g eV, %.4f u, %g q. %.4g points: %s mm",
           _G.ions_flown, _G.ion_energy[group_index], group_mass[group_index], group_charge[group_index], points_used, _G.FLY2_SOURCE_DISTR)
           --_G.ions_flown, _G.ion_energy[group_index], ion_mass, ion_charge, points_used, _G.FLY2_SOURCE_DISTR)
    else
      flown_message = string.format("Flew %d ions @ %g eV, %.4f u, %g q. %.4g points: pattern %d, spacing %.3g mm",
           _G.ions_flown, _G.ion_energy[group_index], group_mass[group_index], group_charge[group_index], points_used, source_point_pattern, _source_point_spacing)
           --_G.ions_flown, _G.ion_energy[group_index], ion_mass, ion_charge, points_used, source_point_pattern, _source_point_spacing)
    end

    local overall_stats
    if math.floor(points_used) == points_used or _isotropic >= 3 then -- Integer number of points (no point missing any directions), or using _isotropic>=3 where the number of points per direction is expected to vary anyway
      print(flown_message)
      
      -- Calculate sum of variances (equal weight for each direction regardless of its hit count)
      y_Q_sum = 0
      t_Q_sum = 0
      y_rel_variance_sum = 0
      y_rel_weight_sum = 0
      t_rel_variance_sum = 0
      t_rel_std_count = 0
      y2_std_sum = 0
      y2_rel_variance_sum = 0
      
      -- Print results for each direction's ion "bunch"
      for dir = 1, _number_of_directions do
        if _isotropic == 2 then
          -- List the directions from +90 to -90 degrees, which visually agrees with the simulation (highest p_Y component on top).
          -- Typically needs odd _number_of_directions (e.g. 5, 9 or 17) because +90 and -90 are not merged.
 
          -- Distinct rays for p_Y>0 and p_Y<0 (to be able to test y-offset sensitivity easily):
          elevation = 90 - (180/(_number_of_directions-1)) * (dir - 1);
          -- Alternative, only at p_Y>0 (to not waste simulation effort in symmetric case):
          --elevation = 90 - (90/(_number_of_directions-1)) * (dir - 1);

        elseif _isotropic >= 3 then
          -- The elevation used has been stored in an array
          elevation = direction_elevation[dir];
        else -- All other methods list directions from 0 to 359 degrees
          elevation = (360/_number_of_directions) * (dir - 1) --[degrees]
        end
        
        abssin_elevation = abs(sin(elevation*DEGREE))
        if det_count[group_index][dir] > 0 then
          ---- Prepare data for the overall statistics of y and t errors:

          -- At least RES_R_ABS [mm] std.deviation (approx. detector resolution) for each bunch
          local y_variance_limited = max(RES_R_ABS*RES_R_ABS, det_y_Q[group_index][dir]/det_count[group_index][dir])
          y_Q_sum = y_Q_sum + y_variance_limited -- [(mm)^2]
          -- At least RES_TOF_ABS [ns] std.deviation (still less than approx. detector resolution 1 ns)
          t_Q_sum = t_Q_sum + max(RES_TOF_ABS*RES_TOF_ABS, det_t_Q[group_index][dir]/det_count[group_index][dir])  -- [(ns)^2]

          -- calculate total relative error as sqrt(sum( rel.err1^2 + rel.err2^2 + ...))
          -- where rel.err.i^2 = (std.i / mean.i)^2 = variance.i/(mean.i*mean.i).
          -- to handle cases where y^2 is near zero (which is OK!), use y=MIN_R_REL (rather large to not get enormous relative error)
          
          -- To avoid dividing by zero or a very small radius, use a (small) threshold MIN_R_REL.
          -- Then the relative error found is scaled by approximately sin(|elevation|),
          -- more precisely: in the RMS (variance sum) the weight for (relative error)^2
          -- is MIN_R_WEIGHT + sin^2(elevation).
          -- Thus the error is (almost) ignored for angles 0 and 180 (where it does not affect energy resolution, only angular resolution).
          local y_limited --[mm] y or a higher value (limit) to divide with to get relative radial "error" (std.y / y)
          y_limited = max(MIN_R_REL, abs(det_y_mean[group_index][dir]))
          local y_rel_weight
          y_rel_weight = MIN_R_WEIGHT + abssin_elevation*abssin_elevation
          y_rel_variance_sum = y_rel_variance_sum + y_rel_weight * max(RES_R_REL*RES_R_REL, det_y_Q[group_index][dir]/det_count[group_index][dir])/(y_limited*y_limited)
          y_rel_weight_sum = y_rel_weight_sum + y_rel_weight

          y2_std_sum = y2_std_sum + sqrt(det_y2_Q[group_index][dir]/det_count[group_index][dir])
          local y2_limited 
          y2_limited = max(MIN_R_REL*MIN_R_REL, abs(det_y2_mean[group_index][dir]))
          y2_rel_variance_sum = y2_rel_variance_sum + y_rel_weight * det_y2_Q[group_index][dir]/det_count[group_index][dir]/(y2_limited*y2_limited)

          dir_front = _number_of_directions/2 + 2 - dir -- same radial velocity component but opposite v_z sign
          -- (need to have dir_front < dir so that all dir_front trajectories are completed)
          t_variance_mean = -1; t_from_center = -1; -- define in case the if-statment here does not run
          if _isotropic ~= 2 and _number_of_directions % 2 == 0 and dir > _number_of_directions/4+1 and 
              iif(dir_front>=1, det_count[group_index][dir_front], 0) > 0 and (dir-dir_front)*(360/_number_of_directions) > 11.25 then
            --DEBUG printf("Comparing %d and %d", elevation, (360/_number_of_directions)*(dir_front-1)) --DEBUG
            t_from_center = max(abs(det_t_mean[group_index][dir_front] - det_t_mean[group_index][dir])/2, 1e-3) --[us] at least 1 ns (just to avoid division by zero)
            -- Use at least RES_TOF_REL [us] standard deviation, otherwise the RMS-average of the two directions' std.:
            t_variance_mean = max(RES_TOF_REL*RES_TOF_REL, ((det_t_Q[group_index][dir]/det_count[group_index][dir]) + (det_t_Q[group_index][dir_front]/det_count[group_index][dir_front]))/2) -- [(us)^2]
            t_rel_variance_sum = t_rel_variance_sum + t_variance_mean/(t_from_center*t_from_center)
            t_rel_std_count = t_rel_std_count + 1
            t_rel_str = string.format("%5.2f%%", 100*sqrt(t_variance_mean)/t_from_center )
          elseif _isotropic == 1 and _number_of_directions == 1 then
            -- A single "direction" group, for isotropic initial angles. This can be used for mass resolution estimate directly from Lua.
            -- For 100u, the mass resolution estimate t/2/std(t) is about 4.1 times as large as
            --  (t_101u-t_100u) / mean([width_101u width_100u]) where width is measured as 5-95% interquantile range.
            t_rel_str = string.format("%4.1fu", 
                det_t_mean[group_index][dir] / 2 / sqrt(det_t_Q[group_index][dir]/det_count[group_index][dir]) )
          else            
            t_rel_str = "" -- no relative TOF resolution calculated for this direction (or _number_of_directions is odd which gives NaN for time-resolution)
          end
          
          ---- Show statistics for this direction
          local fmt
          if _isotropic >= 2 then -- 2D style results (_isotropic = 2, 3 or 4)
            -- (The reduced number of digits is only to suit GUI usage when looking at y^2 column too, could be increased again.)
            fmt = "%3d d&: %2d ions; y %7.3f std %6.3f, range %5.2f mm; t%5.2f std%4.2f ns%s; y2 std%6.2f mm^2 %6.3f%%"
          else -- 3D style results
            fmt = "%3d d.: %2d ions; y %7.3f std %6.3f, range %6.3f mm; t%9.3f std%6.3f ns %s; y2 std%6.2f mm^2 %6.3f%%"
          end
          printf(fmt, math.floor(elevation + 0.5), det_count[group_index][dir],
                det_y_mean[group_index][dir],
                sqrt(det_y_Q[group_index][dir]/det_count[group_index][dir]),
                (det_y_max[group_index][dir] - det_y_min[group_index][dir]),
                det_t_mean[group_index][dir]*1000,
                sqrt(det_t_Q[group_index][dir]/det_count[group_index][dir])*1000,
                t_rel_str, 
                sqrt(det_y2_Q[group_index][dir]/det_count[group_index][dir]), --(note: mean det_y2_mean[group_index][dir] is very near det_y_mean[group_index][dir]^2)
                sqrt(det_y2_Q[group_index][dir]/det_count[group_index][dir])/y2_limited*100 --alternative: sqrt(det_y2_Q[group_index][dir]/det_count[group_index][dir]))/(det_y2_mean[group_index][dir]*det_y2_mean[group_index][dir])
                )

        else
          -- No hit for this direction
          if _isotropic >= 2 then -- 2D style results (_isotropic = 2, 3 or 4)
            printf("%3d d&: %2d ions; y  0      std      0, range  0    mm; t   -0 std  0  ns; y2 std   0    mm^2  0%%", elevation, det_count[group_index][dir])
          else -- 3D style results
            printf("%3d d.: %2d ions; y  0      std      0, range  0     mm; t   -0       std  0     ns ; y2 std   0    mm^2  0%%", elevation, det_count[group_index][dir])
          end
        end
      end

      local total_rel_err = sqrt((y_rel_variance_sum/y_rel_weight_sum + t_rel_variance_sum/t_rel_std_count)/2)
      
      overall_stats = string.format(
          "Tot sy %6.3f mm, <sy/y> %6.3f%%. st %6.3f ns, <st/t> %6.3f%%. <1/R> %6.3f%%",
          sqrt(y_Q_sum/_number_of_directions), -- root(mean(variance(y_limited)))
          sqrt(y_rel_variance_sum/y_rel_weight_sum) * 100,
          sqrt(t_Q_sum/_number_of_directions)*1000,
          sqrt(t_rel_variance_sum/t_rel_std_count) * 100,
          total_rel_err * 100) -- this is where old (before 2015) output ended
      overall_stats = string.format("%s\nTot s(y^2) %6.3f mm^2, <s(y^2)/y^2> %7.4f%%", overall_stats,
          y2_std_sum/_number_of_directions, -- note this is mean(std(y^2)), not root(mean(variance(y^2)))
          sqrt(y2_rel_variance_sum/y_rel_weight_sum) * 100 )
  
      --if ion_number == _G.ions_flown_all_groups then
      if is_last then
        -- Only for last ion (in last group)
		--print(string.format("DEBUG: last with %d. Could we use copy %d ", ion_number, copy_of_ion_number))
        
        -- Compute approximative mass resolution for each side
        if B_MCP > 0 and U[A_MCP] ~= nil and U[B_MCP] ~= nil then
          -- Don't try to compute before all potential array instances have been visited (so U is fully populated)
          
          local mass_max_A, mass_max_B = mass_resolution(_G.ion_energy[group_index]);
          local mass_max_A_02, mass_max_B_02 = mass_resolution(0.2);
          local mass_stats = string.format(
              "\nApproximative max mass resolved A %.0f u, B %.0f u at %.1f eV,   A %.0f u, B %.0f u at %.1f eV",
              mass_max_A, mass_max_B, _G.ion_energy[group_index],
              mass_max_A_02, mass_max_B_02, 0.2);
      
          --if ion_mass > 0.5 or write_DLT == 0 or write_DLT >= 3 then
          if group_mass[group_index] > 0.5 or write_DLT == 0 or write_DLT >= 3 then
            -- Not for electrons when writing first or middle part in DLT-file, but otherwise: append mass resolution to file comment
            overall_stats = overall_stats .. mass_stats
          end
        end
      end
      overall_stats, tmp = string.gsub(overall_stats, "-1[.]#IO%%", "NaN%%") -- show NaN% instead of -1.#IO%
      print(overall_stats)
      
    else
      printf("! ERROR: %d ions flown, giving NON-INTEGER %d points. Some point is missing velocity directions.", _G.ions_flown, points_used)
	  -- If a rather small deviation, this could be due to some particles being initialized by fly2-script to points outside simulation volume 
	  -- so that SIMION doesn't even call segment.initialize() for them.
    end
    
    if write_DLT ~= 0 then
      -- Update global (file) counters
      --_G.DLT_hit_count = _G.DLT_hit_count + _G.ions_flown_all_groups -- TODO: use actual number of writes instead (to allow misses)
      _G.DLT_trigger_count = _G.DLT_trigger_count + _G.ions_flown_all_groups
      if DoNotSetPotentials == 0 then -- the usual potential program
        if write_DLT == 1 or write_DLT == 4 then -- For the first (or only) run: show potentials, clear old comment
          _G.DLT_comment = _G.DLT_comment .. string.format(
                              "U_A: %d %d %d, f %d, p %d, MCP %d,\nU_B: %d %d %d MCP %d", 
                              U[A_FIRST], _U_Afree, U[A_LAST], _A_free, _A_plateau, U[A_MCP],
                              U[B_FIRST], _U_Bfree, U[B_LAST], U[B_MCP])
          if _isotropic == 0 then
            _G.DLT_comment = _G.DLT_comment .. "\nVelocity x,y-distribution: only along y-axis"
          elseif _isotropic == -1 then
            _G.DLT_comment = _G.DLT_comment .. "\nVelocity x,y-distribution: constant phi-density"
          elseif _isotropic == 1 then
            _G.DLT_comment = _G.DLT_comment .. "\nVelocity x,y-distribution: isotropic solid angle density"
          elseif _isotropic == 2 then
            _G.DLT_comment = _G.DLT_comment .. "\nVelocity x,y-distribution: VMI with phi=0 and double weight for theta=90"
          elseif _isotropic == 3 then
            _G.DLT_comment = _G.DLT_comment .. "\nVelocity x,y-distribution: VMI approximately isotropic, including theta=90 (p_Z=0)"
          elseif _isotropic == 4 then
            _G.DLT_comment = _G.DLT_comment .. "\nVelocity x,y-distribution: VMI approximately isotropic, avoiding theta=90 (p_Z=0)"
          else
            _G.DLT_comment = _G.DLT_comment .. "\nVelocity x,y-distribution: " .. _isotropic
          end
        end
      else
        _G.DLT_comment = _G.DLT_comment .. "Customized potentials" --TODO: show all values in this case
      end
      _G.DLT_comment = _G.DLT_comment .. "\n" .. flown_message .. "\n" .. overall_stats:gsub("[<>]","")
    end
    if write_DLT >= 3 then
      -- Finish DLT file & footer
      local foot = "\0\0\0\0\254\254\254\254" --Data end marker
      -- # Counters
      foot = foot .. "\0" .. u24_as_string(_G.DLT_hit_count) --events (groups)
      foot = foot .. "\0" .. u24_as_string(_G.DLT_hit_count) -- accepted hits
      foot = foot .. "\0" .. u24_as_string(_G.DLT_trigger_count) -- start triggings
      foot = foot .. "\0" .. u24_as_string(0) --lone start TODO misses?
      foot = foot .. "\0" .. u24_as_string(0) --discarded events(groups)
      foot = foot .. "\0" .. u24_as_string(0) --discarded complete hits
      foot = foot .. os.date("%Y-%m-%dT%H:%M:%S.000+00:00:00") -- time zone and milliseconds was not supported in this Lua
      foot = foot .. "\0" -- Property list terminator and Comment string terminator, in backward direction
      if _G.DLT_hit_count ~= _G.DLT_trigger_count then
        -- Highlight the fact that some particles didn't hit the detector
        _G.DLT_comment = _G.DLT_comment:gsub("\n", "\nSome particles missed the detector.\n", 1)
      end
      foot = foot .. _G.DLT_comment .. "\0" --Comment String terminator

      output_DLT_file:write(foot)
      output_DLT_file:write("\0" .. u24_as_string(4 + foot:len()))
      output_DLT_file:close()
    elseif output_DLT_file ~= nil then
      output_DLT_file:close()
    end
end

function mass_resolution(energy)
  -- Spectrometer lengths (A- and B-side)
  -- The extraction distance, from source to first mesh
  local D_A = (Z_EXTRACTOR_A - source_point_z) * 0.001; --[m]
  local D_B = (source_point_z - Z_EXTRACTOR_B) * 0.001; --[m]

  local Lb_A = (Z_END_A - Z_EXTRACTOR_A) * 0.001; --[m] the bend region. (drawing says 140 to 141, the latter agrees better with SIMION field)
  local Lf_A = (detectorA_x-Z_END_A) * 0.001;     --[m] the final-acceleration region
  local Lb_B = (Z_EXTRACTOR_B - Z_END_B) * 0.001;
  local Lf_B = (Z_END_B-detectorB_x) * 0.001;

  local U_needle = ((source_point_z-Z_EXTRACTOR_B)*U[A_FIRST] + (Z_EXTRACTOR_A-source_point_z)*U[B_FIRST])/(Z_EXTRACTOR_A-Z_EXTRACTOR_B)
  local V_ext_A = U_needle - U[A_FIRST]
  local V_ext_B = U_needle - U[B_FIRST]

  -- The detector's time resolution (and z width of source) causes some blur,
  -- the detection should be less than 1 ns. I haven't properly estimated what
  -- blur a reasonable source width would give.
  local TOF_FWHM = max(RES_TOF_ABS, RES_TOF_REL) * 1e-6; --[s] FWHM, i.e. more than twice the std of a Gaussian


  local Delta_m = 1; --[u] testing adjacent mass numbers (could change to test like 12 and 16 (CO) to find max energy separable there if more interesting)
  -- The charge is always assumed to be +-1

  -- A-side:
  local kappa = math.abs(energy / V_ext_A); --[eV]/[eV] = dimensionless
  local b = (U[A_FIRST]-U[A_LAST]) / V_ext_A -- bend parameter
  local f = (U[A_LAST]-U[A_MCP]) / V_ext_A -- final acceleration parameter
  
  -- Compared to notes on paper from 2011-04-11, the one (1) has been
  -- subtracted (moved from RHS to LHS):
  local LHS = (    1 /(math.sqrt(kappa+1    ) + math.sqrt(kappa    ))
          + Lb_A/D_A /(math.sqrt(kappa+1+b  ) + math.sqrt(kappa+1  ))
          + Lf_A/D_A /(math.sqrt(kappa+1+b+f) + math.sqrt(kappa+1+b))
        ) / math.sqrt(kappa);
-- Quick, approximative version, error at most 1u for t_FWHM=1ns,
  -- V_ext<1000V, m <= 50, energy=1eV (alt. divide by 4.06 for energy=0eV).
  -- Will be used as starting guess for fzero also in range outside these
  -- parameter ranges.
  -- RHS_approx = (m/Delta_m) *(2.02)*2.03 = m * 4.08/Delta_m;
  -- LHS = RHS_approx ==> 
  local mass_max_A = LHS * Delta_m / 4.1;
  
  
  -- B-side:
  kappa = math.abs(energy / V_ext_B); --[eV]/[eV] = dimensionless
  b = (U[B_FIRST]-U[B_LAST]) / V_ext_B -- bend parameter
  f = (U[B_LAST]-U[B_MCP]) / V_ext_B -- final acceleration parameter
  LHS = (      1    /(math.sqrt(kappa+1    ) + math.sqrt(kappa    ))
         + Lb_B/D_B /(math.sqrt(kappa+1+b  ) + math.sqrt(kappa+1  ))
         + Lf_B/D_B /(math.sqrt(kappa+1+b+f) + math.sqrt(kappa+1+b))
        ) / math.sqrt(kappa);
  local mass_max_B = LHS * Delta_m / 4.1;
  
  return mass_max_A, mass_max_B
end
