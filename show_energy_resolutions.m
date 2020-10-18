% This was a script before (using externally defined variables in optimize_potentials_loop) but changed to function now.
%
% NOTE: expecting flys with a column per particle (mass, energy) and one or multiple rows (potentials), 
% which is is the tranpose of a single-potential result from readsim3 (also transpose of what show_mass_resolution.m uses).
%
function show_energy_resolutions(flys, save_figures, batch_name, quantity)

% Show transverse (2D) y and y^2 resolution (prop. to momentum and kinetic energy resoluton)
% as function of the energy.

% Currently assuming there is only one type of particle (mass, charge, spectrometer-side) in the dataset.

% For each row in flys (same lens settings for all columns)
potentials = transpose(collect_field(flys(:,1), 'potentials', true));
V_repeller = potentials(:, flys(1,1).electrode_index.B_FIRST);
V_extractor = potentials(:,flys(1,1).electrode_index.A_FIRST);
extractor_repeller_ratio = collect_field({flys(:,1).lens_parameters_A}, 'extractor_repeller_ratio');
% extractor_repeller_ratio = V_extractor ./ V_repeller;

if flys(1,1).electrode_index.A_FREEGROUP_COUNT >= 1
  tweak_ratio = potentials(:, flys(1,1).electrode_index.A_FIRST+1) ./ V_repeller;
  tweak_name = 'Tweak ratio';
  tweak_name_subscript = 'tw';
elseif flys(1,1).electrode_index.B_FIRST == 14 && flys(1,1).electrode_index.A_FIRST == 15 && length(potentials) >= 16 ...
      && length(flys(1,1).workbench_name) >= 6 && strcmp(flys(1,1).workbench_name(1:6), 'both8e')
  % For 'both8e', 'both8eu', 'both8eu_B2': although tweak ratio
  tweak_ratio = potentials(:, flys(1,1).electrode_index.A_FIRST+1) ./ potentials(:, flys(1,1).electrode_index.A_FIRST);
  tweak_name = 'Divider 1 = V_{16}/ V_{15}  ';
  tweak_name_subscript = 'Div.1';
else
%   tweak_ratio = []; % original
  tweak_ratio = [NaN]; % 2015-06-23 to try to plot a single flys
  tweak_name = 'not used';
  tweak_name_subscript = '?';
end

which = true(size(flys,1),1); % all lens settings
% which = abs(tweak_ratio) < 0.005;
% which = tweak_ratio > 0.005;
% %which = abs(extractor_repeller_ratio-0.65) < 0.045; %the range around 0.65 seems optimal, when averaging or optimizing other parameters (for 352mm and maybe 202mm)
% %which = abs(tweak_ratio) < 0.005 & abs(extractor_repeller_ratio-0.65) < 0.045; % essentially same range optimal (and as good) for tweak_ratio=0 too. Thus the free electrode is OK to keep grounded for now.
% which = abs(tweak_ratio) < 0.005 & abs(extractor_repeller_ratio-0.61) < 0.035 & V_repeller > -3001; %(for 342mm % 302: optimized ranges to retain max or near max, but maximize average) Selects 12 settings.
% which = abs(tweak_ratio) < 0.005 & abs(extractor_repeller_ratio-0.71) < 0.015; %(for 342mm_iontwak)
% which = abs(V_repeller - -650) <= 1; % a specific repeller potential
% which = abs(V_repeller - -1350) < 1 & abs(extractor_repeller_ratio-0.685)<0.002;
% which = abs(V_repeller - -1500) < 1 & abs(tweak_ratio-0.3) < 0.005;
% which = abs(extractor_repeller_ratio-0.66)<0.002;

% To be able to skip ion and look at electron-VMI only. (NOTE: this modifies the workspace flys structure!)
flys = flys(:, collect_field({flys(1,:).source}, 'mass') < 1);  % only electrons (at all energies)

energies_used = unique(transpose(collect_field({flys(1,:).source}, 'energy')));

% -------------
figure(6); clf;
%set(gcf,'Position',[ 740   191   833   775]);
set(gcf,'Position',[ 440   35   833   775]);
% subplot(4,2,1:4);
subplot(4,2,[1 3]);
indices = find(which);

% Called to show figure for a single lens setting, rather than with a large set of potentials?
single_setting_figure_even_if_misses = (size(flys,1) == 1 && length(indices) == 1 && strcmp(batch_name, 'just evaluate initial simulation'));

all_pow_y = zeros(size(flys));
all_pow_y2 = zeros(size(flys));
angles = flys(1,1).rays(:,strcmp(flys(1,1).ray_columns, 'theta'))';
angular_pow_y2 = NaN(size(flys(1,1).rays,1), size(flys,2), length(indices)); % 3D matrix: (angle, energy, index)

counter = 0;
for i = indices'
  fly = flys(i,:);
  
  misses = [fly.misses];

  if ~any(misses) ... % (the original: drop a lens setting that gave some miss)
      || (single_setting_figure_even_if_misses && sum(misses==0) >= 2) % or at least two energies are without misses
    summary = collect_field(fly,'overall');
    
    % relative y^2 and y
    all_pow_y(i,:)  = 1./[summary.y_rel];
    all_pow_y2(i,:) = 1./[summary.y2_rel];
    if single_setting_figure_even_if_misses
      all_pow_y(i, misses>0) = 1; % leave a non-NaN value to keep energy axis to full range for easier comparison to flys without misses
      all_pow_y2(i, misses>0) = NaN;
    end
    h = plot(energies_used, all_pow_y(i,:), '--r', energies_used, all_pow_y2(i,:), ':b');
    hold on;
    counter = counter + 1;

    % Another version of y^2 power, averaging all energies but distinguishing the angles
    rays = [fly.rays]; % give a N_DIRECTIONS-by-(10*size(flys,2)) matrix.
    if all(size(fly(1).source.spacing) == [3 1]) % this (not transposed) identifies the 2015-04-02 version where std/y2 ratios were written incorrectly by Lua
      % variant to correct the error
      angular_pow_y2(:,:,counter) = 1 ./ sqrt( rays(:, find(strcmp(fly(1).ray_columns,'y2_std_rel')):length(fly(1).rays):end) );
    else % normal
      angular_pow_y2(:,:,counter) = 1 ./ rays(:, find(strcmp(fly(1).ray_columns,'y2_std_rel')):length(fly(1).rays):end);
    end
    if single_setting_figure_even_if_misses
      angular_pow_y2(:, misses>0, counter) = NaN;
    end
    
%     subplot(3,1,2);
%     plot(energies_used, [summary.y2_std_mean], '-b')
%     hold on;
%     
%     subplot(3,1,3);
%     plot(energies_used, [summary.y_std], '-r')
%     hold on;
  else
    % Remove this simulation from selection, because not all energies were fully detected
    which(i) = false;
  end
end
ylabel('Relative resolving power');
xlabel('Kinetic energy / eV')
set(gca,'Position', [0.13 0.58 0.775 0.4])
% set(gca,'Position', [0.08 0.6 0.4 0.4])

number_of_preselected = length(indices);
indices = find(which);

angular_pow_y2(:,:,counter+1:end) = []; % remove extra preallocated entries, if any
%until 2015-12-18 11:59: mean_pow_y  = mean(all_pow_y(which,:), 2);
%until 2015-12-18 11:59: mean_pow_y2 = mean(all_pow_y2(which,:),2);
% 2-norm average the 1/power among energies, to put more weight on worst energy, striving to get all OK rather than one very good and others worse:
% mean_pow_y  = 1./sqrt(mean(1./all_pow_y(which,:).^2, 2));
% mean_pow_y2 = 1./sqrt(mean(1./all_pow_y2(which,:).^2,2));
% 3-norm average the 1/power among energies, to put even more weight on worst energy, striving to get all OK rather than one very good and others worse:
mean_pow_y  = 1./mean(1./all_pow_y(which,:).^3, 2).^(1/3);
mean_pow_y2 = 1./mean(1./all_pow_y2(which,:).^3,2).^(1/3);


weights = abs(sin(angles*pi/180)); weights = weights / sum(weights);
% mean_pow_y2_sin_weighted = (weights * squeeze(mean(angular_pow_y2,2)))'; % averaging across energies, then weighted average across angles ==> scalar per lens setting
all_pow_y2_sin_weighted = squeeze(sum(repmat(weights',[1 size(angular_pow_y2,2) size(angular_pow_y2,3)]) .* angular_pow_y2,1))'; % weighted average across angles ==> scalar per (lens setting, energy) combination
if size(all_pow_y2_sin_weighted,2) == 1
  % When only a single fly is selected, the squeeze turns it into a column vector while it should be a row per lens setting
  all_pow_y2_sin_weighted = transpose(all_pow_y2_sin_weighted);
end

mean_pow_y2_angle_filtered = squeeze(mean(angular_pow_y2(abs(angles)>25,:,:),1))'; angular_filter = '|theta|>25'; % only large angles (more important for VMI and easier to focus)
% mean_pow_y2_angle_filtered = squeeze(mean(angular_pow_y2(abs(angles)<=25,:,:),1))'; angular_filter = '|theta|<=25'; % only small angles

if isempty(all_pow_y2_sin_weighted)
  h_filtered = [NaN];
  h_sin = [NaN];
else
  h_filtered = plot(energies_used, transpose(mean_pow_y2_angle_filtered), '-g');
  h_sin = plot(energies_used, transpose(all_pow_y2_sin_weighted), '-b');
end

%legend('y_{limited} / std(y)', 'y_{limited}^2 / std(y^2)', 'Location','best');
set(gca,'yscale','log')
if isempty(indices)
  text(0.1, 100, 'Only misses');
  text(0.1, 0.7, 'Only misses');
else
  legend([h; h_filtered(1); h_sin(1)], 'y_{limited} / std(y)','y_{limited}^2 / std(y^2)',['-''''- avg, ' angular_filter],'-''''- |sin\theta|-weighted avg.', 'Location','EO');
end

% axis tight
% yl = ylim(); yl(1) = max(yl(1),4); ylim(yl);
% ylim([4 150])
% ylim([8 300]) % used mostly until 2015-12-15
%ylim([8 300] * 16/8) % used from 2015-12-15
ylim([6 500]) % 2020-02-26. For interst in low energy electrons and prioritizing mass side, even moderate resolvig powers worth seeing

% subplot(4,2,[2 4]);
% hold on
% h = plot(energies_used, mean_pow_y2_angle_filtered, '-g', energies_used, all_pow_y2_sin_weighted, ':b');
% xlabel('Kinetic energy')
% ylabel('y_{limited}^2 / std(y^2)');
% set(gca,'Position', [0.58 0.6 0.4 0.3768])
% legend(h([1 end]), [angular_filter ' avg.'], '|sin\theta|-weighted avg.', 'Location','best');

subplot(4,2,5:6);

%until 2015-12-18 11:59: mean_pow_y2_sin_weighted = mean(all_pow_y2_sin_weighted,2);
%until 2016-01-06 15:50: mean_pow_y2_sin_weighted = 1/sqrt(mean(1./all_pow_y2_sin_weighted.^2,2)); % 2-norm average the 1/power among energies, to put more weight on worst energy, striving to get all OK rather than one very good and others worse
mean_pow_y2_sin_weighted = 1/mean(1./all_pow_y2_sin_weighted.^3,2).^(1/3); % 3-norm average the 1/power among energies, to put even more weight on worst energy, striving to get all OK rather than one very good and others worse

if isempty(mean_pow_y)
  badness_mix_as_in_optimizer = [];
else
  % hist(mean_pow_y2, 20) % old before 2015-09-04
  %until 2015-12-18 11:59: badness_mix_as_in_optimizer = (1./mean_pow_y + 1./mean_pow_y2 + 2./mean_pow_y2_sin_weighted) / 4; % Average (already energy-averaged) quantities with double weight on y2_rel_sin
  %until 2015-12-30 : badness_mix_as_in_optimizer = (2./mean_pow_y + 1./mean_pow_y2 + 2./mean_pow_y2_sin_weighted) / 5; % Average (here already energy-averaged, while optimizer does energy 2-norm average at end) quantities with double weight on y2_rel_sin and y_rel, to typically avoid high-scoring when y2_rel optimized for lowest energy at the cost of other
  badness_mix_as_in_optimizer = (1./mean_pow_y + 1./mean_pow_y2 + 1./mean_pow_y2_sin_weighted) / 3; % Average quantities with equal weight. In recent improvement of Lua code, it was discovered that no y-center coordinate was used (thus falsly good resolution when detector off-centre by 87 mm). Now it seems OK to average with equal weight
end

power_mix = 1 ./ badness_mix_as_in_optimizer; % one row with a combined "resolving power" per evaluated configuration
hist(power_mix, 20);
if length(power_mix) == 1
   xlabel(sprintf('Mixed resolving power from RMS = %.1f', power_mix));
else
   xlabel('Mixed resolving power from RMS'); % new 2015-09-04
end

[max_pow_y2_sin_weighted, best_index_among_indices] = max(mean_pow_y2_sin_weighted);
title(sprintf('{\\bfAverages} among %d selected: V_{rep}=%.0f, \\rho_r=%.3f, \\rho_{%s}=%0.3f\\pm%0.1f\ny^2/std = %.1f (max %.1f), sin\\theta-weighted %.1f (max %.1f); y/std = %1.f (max %.1f)', ...
  length(indices), mean(V_repeller(which)), mean(extractor_repeller_ratio(which)), tweak_name_subscript, mean(tweak_ratio(which)), std(tweak_ratio(which)), ...
  mean(mean_pow_y2), max(mean_pow_y2), mean(mean_pow_y2_sin_weighted,1), max_pow_y2_sin_weighted, ...
  mean(mean_pow_y), max(mean_pow_y) ));
% the maxima are for row-averages, i.e. a lens setting that gave maximal average of resolving power among energies
set(gca,'Position', [0.13 0.3291 0.7750 0.13])

if length(indices) <= 30
  % To show selectd potentials
  disp(sprintf('Printing %d all selected & successful settings:', length(indices)))
  indices_to_show = 1:length(indices);
else
  disp('Only printing best setting:')
  indices_to_show = best_index_among_indices;
end
for ii = indices_to_show
  i = indices(ii);
  disp(sprintf('(%4d) Potentials %5.0f %5.0f %4.0f; Resolv.pow. y %3.0f, y^2 %5.1f, y^2sin %5.1f', ...
    i, V_repeller(i), potentials(i,flys(1,1).electrode_index.A_FIRST), potentials(i,flys(1,1).electrode_index.A_FIRST+1), ...
    mean_pow_y(ii), mean_pow_y2(ii), mean(all_pow_y2_sin_weighted(ii,:),2)));
end
% 202mm: index 353 & 245 & 299, 302mm: 242, 352mm: 296 & 242, 352iontweak: 365


% Scatter plot with ring size showing resolving power
ring_scale = 16 / max(mean_pow_y2);
subplot(4,2,7); cla;
hold on;
for i = 1:length(indices)
  plot(extractor_repeller_ratio(indices(i)), V_repeller(indices(i)), 'o', 'MarkerSize', max([1 mean_pow_y2(i)*ring_scale])); 
end
%xlim([floor(min(extractor_repeller_ratio(indices(i)))/0.1), ceil(max(extractor_repeller_ratio(indices(i)))/0.1)]*0.1);
xlim([min([0.5 min(extractor_repeller_ratio(indices(i)))]), 1])
xlabel('Extractor ratio');
ylabel('V_repeller', 'Interpreter','none');

subplot(4,2,8); cla;
hold on;
for i = 1:length(indices)
  plot(extractor_repeller_ratio(indices(i)), tweak_ratio(indices(i)), 'o', 'MarkerSize', max([1 mean_pow_y2(i)*ring_scale])); % very clearly correlated
%   plot(V_repeller(indices(i)), tweak_ratio(indices(i)), 'o', 'MarkerSize', max([1 mean_pow_y2(i)])); %-- not so clearly correlated
end
xlim([floor(min(extractor_repeller_ratio(indices(i)))/0.1), ceil(max(extractor_repeller_ratio(indices(i)))/0.1)]*0.1);
%ylim([floor(min(tweak_ratio(indices(i)))/0.1), ceil(max(tweak_ratio(indices(i)))/0.1)]*0.1);
ylim([floor(min(tweak_ratio(indices(i)))/0.1)*0.1, max([ceil(max(tweak_ratio(indices(i)))/0.1)*0.1 0.7])]);
xlabel('Extractor ratio');
ylabel(tweak_name);


% To show selected parameters:
% [V_repeller(which)/1000 extractor_repeller_ratio(which) tweak_ratio(which)]

selection_name = '';
if all(diff(V_repeller(which)) == 0)
  selection_name = sprintf('%s rep%.0f', selection_name, mean(V_repeller(which)));
end
if abs(diff(minmax(extractor_repeller_ratio(which)'))) <= 0.021
  selection_name = sprintf('%s re%.3g', selection_name, mean(extractor_repeller_ratio(which)));
end
if abs(diff(minmax(tweak_ratio(which)'))) <= 0.001
  selection_name = sprintf('%s rt%.3g', selection_name, mean(tweak_ratio(which)));
end

if exist('save_figures','var') && save_figures % NOTE: Using variable set by show_simulated_per_energy

  if single_setting_figure_even_if_misses
    % Called from optimze_potentials_loop, for a single flys
    potentials_str = ['-' strrep(strrep(sprintf(' %+.0f', round(flys(1).potentials')), '-',''), '+0','0')]; % skip (1) which is the grounded walls
    potentials_str = strrep(potentials_str, '- 0', '- '); % skip potential (1) which is the grounded walls (even if no + sign before the zero)
    
    if exist('workbench_name','var')
      potentials_str = [workbench_name ' ' potentials_str];
    end
    if isempty(strfind(char(quantity), 'mass')) && isempty(strfind(char(quantity), 'iso'))
      % Add a prefix when using eVMI alone and without the new almost-isotropic distribution
      potentials_str = ['only ' potentials_str];
    end
    
    if length(flys(1).source.spacing) == 3 % 3D Gaussian
      potentials_str = sprintf('%s W%.1f', potentials_str, flys(1).source.spacing(2) * (2*sqrt(2*log(2)))); %[mm]
%       annotation_suffix = [annotation_suffix sprintf('. Source FWHM: flight axis %.1f, optical axis y %.1f, vertical %.1f [mm].', ...
%                 flys(1).source.spacing(1) * (2*sqrt(2*log(2))), flys(1).source.spacing(2) * (2*sqrt(2*log(2))), flys(1).source.spacing(3) * (2*sqrt(2*log(2))) )]; %[mm]
    end
    
    if round(flys(end).source.energy) == 20
      export_fig_multiformat(6, sprintf('VM20 %s', potentials_str), 'PNG');
    elseif round(flys(end).source.energy) <= 7
      export_fig_multiformat(6, sprintf('VM7 %s', potentials_str), 'PNG');
    else % default mostly used for 30 eV
      export_fig_multiformat(6, sprintf('VMI %s', potentials_str), 'PNG');
    end
    
  else % original
    export_fig_multiformat(6,sprintf('resolutions %s%s sel%d', batch_name, selection_name, number_of_preselected), 'PNG');%,'PDF')
  end
end

% Log 2015-04-02
%  Looking at y^2 resolving power for Eppik_style with lengths 202, 302 and 352 mm. (No filtering here.)
%  The max resolving power (averaged across all energies, as well as visual "envelope") increases with length, but not very much.
%  The short length allows lower repeller voltages, and generally requires higher extractor ratio (stronger lens).
%   - Conversely, the long length disfavours low repeller voltages and requires lower extractor ratio (weaker lens).
%  In the unfiltered dataset, the average (of average across energies) is highest in 302 mm, perhaps due to effects of misses when long weak lens, and aberration when short and strong lens?
%  However, if filtering to look only at tweak_ratio=0 ("sel144") then length 202 is best, due to gain in the mid-energies.
% ----

