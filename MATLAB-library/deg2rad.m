% angle_in_radians = deg2rad(angle_in_degrees)
% 
% Convert angle from [radian] unit to [degree] unit.
%
% PARAMETERS
%   angle_in_degrees
% RETURNS
%   angle_in_radians
% SEE ALSO
%   rad2deg rad2principal
function angle_in_radians = deg2rad(angle_in_degrees)

angle_in_radians = angle_in_degrees * pi/180;