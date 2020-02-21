function [euler, qrot] = rodmis2eu(rod,q1)
%RODMIS2EU Convert Rodrigues vector misorientation to euler angles
%%% Inputs:
% -rod: Rodrigues vector of misorientation as 1x3 array
% -q1: Quaternion of orientation corresponding to R=[0,0,0] as 1x4 array
%%% Outputs:
% -euler: Euler angles in Bunge convenction as 1x3 array
% -qrot: Quaternion of orientation corresponding to euler

    x = rod(1);  % rodrigues vector
    y = rod(2);
    z = rod(3);
    w = 1/sqrt(1+x^2+y^2+z^2);          % Calculate first quaternion entry
    q2 = [w w*x w*y w*z];   % Quaternion of grid point (centered around north pole)
    qrot = [q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)-q1(4)*q2(4), ...
        q1(1)*q2(2)+q1(2)*q2(1)+q1(3)*q2(4)-q1(4)*q2(3), ...
        q1(1)*q2(3)+q1(3)*q2(1)-q1(2)*q2(4)+q1(4)*q2(2), ...
        q1(1)*q2(4)+q1(4)*q2(1)+q1(2)*q2(3)-q1(3)*q2(2)];   % quaternion grid rotated to be centered on starting Euler angles
    euler = qu2eu(qrot)*180/pi;     % Euler angles of grid point (in degrees)

end

