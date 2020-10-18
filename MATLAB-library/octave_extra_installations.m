 % To install missing GNU Octave packages
 % The UHH laptops had GNU Octave 5.1.0 (2019) in a Windows GUI version (CLI exists too).
 
% https://octave.org/doc/interpreter/index.html#SEC_Contents
% https://octave.sourceforge.io/docs.php
% https://octave.sourceforge.io/
% https://www.gnu.org/software/octave/missing.html


pkg install -forge octave
pkg install -forge io
pkg install -forge statistics

% At each run:
% Octave setup
addpath('(your path to)\Erik''s MATLAB-library')
pkg load statistics
