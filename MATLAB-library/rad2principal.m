% principal_angle = rad2principal(angle)
% 
% Convert any angle to its equivalent angle in the interval
% -pi to +pi radians.
%
% PARAMETERS
%   angle            (radians)
% RETURNS
%   principal_angle  (radians)
% SEE ALSO
%   deg2principal
function principal_angle = rad2principal(angle)

% Simple, but returns -pi instead of +pi for the switching point
% principal_angle = mod(angle+pi,2*pi) - pi;

% Returns +pi for the switching point
principal_angle = mod(angle,2*pi);
wrap = principal_angle > pi
principal_angle(wrap) = principal_angle(wrap) - 2*pi;
