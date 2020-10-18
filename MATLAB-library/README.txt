The files here are originally from Erik Månsson's "MATLAB-library" directory.

It is intended to be included in the Matlab search path so that they
can be used as functions/subroutines in all other matlab programs/scripts.
To do this in Matlab 2014-2017: 
* In the "Home"-ribbon, click the "Set Path" buttton
* Click "Add Folder..."
* Select this directory (called "MATLAB-library" or something like "tools from Erik").

You can optionally also add the subfolder "export_fig", which has been downloaded
from Oliver Woodford (https://github.com/altmany/export_fig) and is useful
for exporting Matlab figures to PDF/EPS/PNG from the command line (better than the print command in Matlab 2014).
In Matlab 2017 the internal export seems good and it is more difficult to make export_fig work (depending on Ghostscript)
so it is not recommended anymore. The export_fig_multiformat.m is no longer calling export_fig.

Ask Erik for clarifications beyond the documemtation that is usually
available within each file, e.g. by typing "help AICc" or "doc AICc"
for the AICc function.

This directory is to be put on your MATLAB path so you can call a function by name
from a script you run in any other directory. There are no example scripts in the present directory.
