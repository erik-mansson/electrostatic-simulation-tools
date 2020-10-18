% export_fig_multiformat(figure_handle, filename_base, formats_and_options ...)
% Save a figure as one or several image files.
%
% PARAMETERS
% figure_handle  which figure to save: result from gcf or an integer.
%
% filename_base  string without filename suffix, possibly with path.
%
% Further arguments are options that may come in any order:
%   If any of 'PNG', 'EPS', 'PDF', 'EMF' or 'FIG' is given as option, 
%   only that/those format is used. Default is to use PNG & PDF.
%   'bg' followed by a colorspecification sets the background color,
%   default is 'white' (but if only saving to PNG the color is not set).
%
% EXAMPLE
%  export_fig_multiformat(1, 'something', 'PDF', 'PNG')
%  to save Figure 1 as something.pdf and something.png
%
function export_fig_multiformat(varargin)

if nargin < 2
  error('Both figure handle and filename base must be given.')
end
figure_handle = varargin{1};
filename_base = varargin{2};
options = varargin(3:end);

% Replace characters used in shell (glob) patterns. "*" --> "#" and "?" --> "~"
if filename_base(2:3) == ':\'
  % Looks like an absolute Windows path, e.g. C:\... or D:\...
else
  % Replace colon (not valid Windows filename) with semicolon or space 
  filename_base = strrep(filename_base, ':',';');
  % filename_base = strrep(filename_base, ':',' ');
end
[p,f,e] = fileparts(filename_base);
f = [f e]; e = ''; % if filename has multiple periods, we want to treat it all as name (no extension)
f = strrep(strrep(f, '*','#'), '?','~'); % replace some other characters in filename part for convenience (not in path)
f = strrep(f, '  ',' '); % avoid double spaces
% NOTE: / and \ remain, so a directory path can be specified as part of the filename
% Rebuild file path
filename_base = fullfile(p, f);
background_color = get_argument_from_cells(options, 'bg', 'val', 'white');

% If any of 'PNG', 'EPS' or 'PDF' is given as option, only that/those
% format is used.
formatlist_PNG = get_argument_from_cells(options, 'PNG', 'bool');
formatlist_EPS = get_argument_from_cells(options, 'EPS', 'bool');
formatlist_PDF = get_argument_from_cells(options, 'PDF', 'bool');
formatlist_EMF = get_argument_from_cells(options, 'EMF', 'bool');
formatlist_FIG = get_argument_from_cells(options, 'FIG', 'bool');
if not(any([formatlist_PNG formatlist_EPS formatlist_PDF formatlist_EMF formatlist_FIG]))
  % If no format is listed, use PNG & PDF
  formatlist_PNG = true;
  formatlist_EPS = false;
  formatlist_PDF = true;
  formatlist_EMF = false;
  formatlist_FIG = false;
end

set(figure_handle,'PaperPositionMode','auto')
if formatlist_PDF || formatlist_EPS
  % if necesary to get white background in PDF
  set(figure_handle,'InvertHardcopy','on', 'Color',background_color); 
  pause(0.05);
end

if formatlist_PNG
  print(figure_handle, '-dpng', '-r0', sprintf('%s.png', filename_base));
end
if formatlist_EPS
  print(figure_handle, '-depsc2', '-r0', sprintf('%s.eps', filename_base));
end
if formatlist_PDF
  pause(0.05);
  % PDF fails if filename contains '%'
% Was used before (Maybe Matlab 2014, Windows 7) but no longer working:  export_fig(sprintf('%s.pdf', strrep(filename_base,'%','')), '-q101', figure_handle);
% disp('DEBUG: using workaround in export_fig_multiformat since ghostscript just gives blank on Windows 10 laptop')
  % If PDF saving keeps using a too small (portrait A4) papersize by the saveas or GUI, and just blank because some Ghotscript failure on Windows 10 laptop,
  % then I guess I should implement the method from
  %  https://uk.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
  % (or with less space but perhaps dependent on axis tick values then https://uk.mathworks.com/help/matlab/creating_plots/save-figure-with-minimal-white-space.html )  
  % New attempt 2018-06-22
%   fig = figure(figure_handle); or gcf() which in Matlab2017 returns matlab.ui.Figure object
%   fig.PaperPositionMode = 'auto'
%   fig_pos = fig.PaperPosition;
%   fig.PaperSize = [fig_pos(3) fig_pos(4)];
  set(figure_handle,'Units','centimeters');
  pos = get(figure_handle,'Position');
  set(figure_handle,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)],'Renderer','painters')
  % NOTE: The key to getting real vector graphics export is to use renderer 'painters' (instead of opengl or zbuffer)!!
  print(figure_handle, sprintf('%s.pdf', strrep(filename_base,'%','')), '-dpdf','-r0');
  set(figure_handle,'Units','pixels');

end
if formatlist_EMF
  print(figure_handle, '-dmeta', '-r0', sprintf('%s.emf', filename_base));
end

if formatlist_FIG % save figure in Matlab's .fig-format
  saveas(figure_handle, sprintf('%s.fig', filename_base), 'fig');
end