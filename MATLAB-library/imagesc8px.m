% h = imagesc8px(x_minmax, y_minmax, img)
% h = imagesc8px(x_minmax, y_minmax, img, bin_width)
% h = imagesc8px(x_minmax, y_minmax, img, bin_width_x, bin_width_y)
%
% Upsample a matrix by 8x8 to make the resulting image
% survive JPEG DCT compression when Inkscape exports to EPS.
%
function h = imagesc8px(x_minmax, y_minmax, img, bin_width, bin_width_y)
if length(x_minmax) > 2
  x_minmax = [min(x_minmax) max(x_minmax)];
end
if length(y_minmax) > 2
  y_minmax = [min(y_minmax) max(y_minmax)];
end

if nargin < 4
  % If no bin width given, compute using axes and number of bins
  bin_width   = abs(diff(x_minmax([1 end]))) / (size(img,2)-1);
  bin_width_y = abs(diff(y_minmax([1 end]))) / (size(img,1)-1);
  
elseif nargin < 5
  bin_width_y = bin_width;
end

im = imresize(img, 8, 'nearest');

% offset from center in the new 8 pixels: [-7 -5 -3 -1 1 3 5 7]/16 px
h = imagesc(x_minmax+[-7 7]*bin_width/16, y_minmax+[-7 7]*bin_width_y/16, im);

set(gca, 'YDir', 'normal'); % let y increase upwards