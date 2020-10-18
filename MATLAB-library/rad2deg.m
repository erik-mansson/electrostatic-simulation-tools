% angle_in_degrees = rad2deg(angle_in_radians)
% 
% Convert angle from [degree] unit to [radian] unit.
%
% PARAMETERS
%   angle_in_radians
% RETURNS
%   angle_in_degrees
% SEE ALSO
%   deg2rad rad2principal
function angle_in_degrees = rad2deg(angle_in_radians)

angle_in_degrees = angle_in_radians * 180/pi;