% Plot (x,y) pairs with errobars along y.
% Created 2015-Aug-27 by Erik Månsson (function-ifying older code snippets)
%
% Styles are currently not supported for the 'serif' type.
%
% PARAMETERS
%  x      abscissa array
%  y      ordinate array
%  l      errorbar (half)length to add or subtract from y to get its endpoints
%  style  linestyle string, the colour should be the first character (otherwise default errorbar_color is gray)
%  errorbar_color (Default: from style or [0.5 0.5 0.5] gray)
%  errorbar_type  (Default: 'area only')
%                'area'       shaded area connecting error"bars" of nearby data points, with edge line
%                'area only'  shaded area connecting error"bars" of nearby data points, without edge line
%                'sans'       simple vertical line
%                'serif'      vertical line with horizontal marks at the endpoints
%  errorbar_alpha (Default: 0.2) opacity. 0 is transparent, 1 is fully opaque
% RETURNS
%  h            Handle to the main curve
%  h_errorbars  Handle to the errorbar(s)
function [h, h_errorbars] = plot_with_errorbars(x,y,l,style, errorbar_color, errorbar_type, errorbar_alpha)
if ischar(l)
  error('Third argument needs to be an array of errorbar (half)lengths. Linestyle is fourth argument.');
end

if nargin < 4
  style = 'k-';
  % IMPROVEMENT: is there a property in Figure or Axes to get list of default colors and cycle through them?
end
if ischar(style)
  if isempty(style)
    style_color = ''; % use default?
    style_color_array = [];
  elseif isempty(strfind('bgrcmykw', style(1)))
    % First character is not a color character.
    style_color = ''; % use default?
%     style_color = 'k'; % use black
    style_color_array = [];
  else
    style_color = style(1);
    style_color_array = style_color; % string is OK also where colour-array accepted
  end
  if strcmp(style_color, 'g') % change green from bright to darker, for better contrast versus white background
    style_color_array = [0 0.8 0]; 
  end
  style_marker = style(2:end);
  style_marker = strrep(strrep(strrep(style_marker, '-.',''), '-',''), ':','');
  style_line   = strrep(style(2:end), style_marker, '');
elseif isnumeric(style) && length(style) == 3 % only colour array
  style_color = style;
  style_marker = '';
  style_line   = '-';
else
  error('Invalid style.');
end
if isempty(style_line)
  style_line = 'none';
end

if nargin < 5 || isempty(errorbar_color)
  errorbar_color = style_color;
  if isempty(strfind('bgrcmykw', errorbar_color))
    % Not a color character. Use a default gray
    errorbar_color = [0.5 0.5 0.5];
  end
  if strcmp(errorbar_color, 'g') % change green from bright to darker, for better contrast versus white background
    errorbar_color = [0 0.8 0]; 
  end
end
if nargin < 6
%   errorbar_type = 'area'; % was default until 2016-05-11
  errorbar_type = 'area only';
end
if nargin < 7
  errorbar_alpha = 0.2;
end

% Convert arrays to column arrays
x = x(:);
y = y(:);
l = l(:);
if isempty(l) % no error bar data
  l = NaN(size(y));
end
if length(x) ~= length(y)
  error('Sizes of x and y arrays differ.');
end
if length(l) ~= length(y)
  error('The array with errorbars lengths and the y array have different size.');
end

switch errorbar_type
  case 'serif' % with horizontal markers at ends of errorbars
    h_errorbars = [];
    h = errorbar(x, y, l, l);
    set(h, 'LineStyle', style_line); 
    set(h, 'Color', style_color_array); 
    set(h, 'Marker', style_marker); 
    
  case 'sans' % simple vertical lines, without horizontal markers at ends of errorbars
    h_errorbars = plot([x x]', [y-l, y+l]', style);%, 'Alpha',errorbar_alpha);
    set(h_errorbars, 'Marker','none', 'Color', errorbar_color);
    if isempty(style_line) || strcmp(style_line, 'none')
      % Errobars always need some line. Use '-' as default
      set(h_errorbars, 'LineStyle','-');
    end
      
    hold on;
    % The main line
    h = plot(x, y, style, 'LineWidth',1.5);
    set(h, 'Color', style_color_array); 
  
  case {'area', 'area only'} % filled area for errorbars

    % A patch may not have any NaN coordinate. Trim and split to plot as much as possible in case there are NaNs present!
    ok = ~isnan(x+y+l);
    starts = ok; starts([false; ok(1:end-1)]) = false; starts = find(starts); % a new block starts where a true value is preceded by a false
    if strcmp(errorbar_type, 'area only')
      edge_color = 'none';
    else
      edge_color = errorbar_color;
    end
    
    for start = starts(:)'
      which = start:(start-2+find([~ok(start:end); true],1,'first'));
      ex = [x(which); flipud(x(which))]';
      ey = [y(which)-l(which); flipud(y(which)+l(which))]';
      h_errorbars = patch(ex, ey, zeros(1,length(ex)), 'FaceColor',errorbar_color, 'EdgeColor',edge_color, ...
                          'LineWidth',0.5, 'FaceAlpha',errorbar_alpha, 'EdgeAlpha',errorbar_alpha);
      hold on;
    end
    % The main line
    h = plot(x, y, style, 'LineWidth',1.5);
    set(h, 'Color', style_color_array); 
    
  otherwise
    error('Unsupported errorbar_type "%s"', errorbar_type);
end