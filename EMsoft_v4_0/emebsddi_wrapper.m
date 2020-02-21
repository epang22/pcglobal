function dp = emebsddi_wrapper(x_all, L, xpc, ypc, q1, Lfactor, xpcfactor, ypcfactor, eulerfactor, options)
%EMEBSDDI_WRAPPER Calculate dot products from scaled optimization inputs
% Original: 7/24/19 (Edward Pang, MIT)
% Change log:
% 12/17/19 ELP: specify optimal nnk, numdictsingle, numexptsingle here

    % Compute euler angles from Rodrigues vector/eulerfactor
    x = x_all(4)*eulerfactor;    % Rodrigues vector components (grid centered about origin)
    y = x_all(5)*eulerfactor;
    z = x_all(6)*eulerfactor;
    w = 1/sqrt(1+x^2+y^2+z^2);          % Calculate first quaternion entry
    q2 = [w w*x w*y w*z];   % Quaternion of grid point (centered around north pole)
    qrot = [q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)-q1(4)*q2(4), ...
        q1(1)*q2(2)+q1(2)*q2(1)+q1(3)*q2(4)-q1(4)*q2(3), ...
        q1(1)*q2(3)+q1(3)*q2(1)-q1(2)*q2(4)+q1(4)*q2(2), ...
        q1(1)*q2(4)+q1(4)*q2(1)+q1(2)*q2(3)-q1(3)*q2(2)];   % quaternion grid rotated to be centered on starting Euler angles
    euler = qu2eu(qrot)*180/pi;     % Euler angles of grid point (in degrees)
    
    % Send to emebsddi_dp.m to compute dot product
    L2 = x_all(1)*Lfactor + L;
    xpc2 = x_all(2)*xpcfactor + xpc;
    ypc2 = x_all(3)*ypcfactor + ypc;
    
    options.nnk = 1;
    options.numdictsingle = 32;   % number of dictionary files arranged in column for dp on GPU (multiples of 16 perform better)
    options.numexptsingle = 32;   % number of experiment files "
    dp = emebsddi_dp(L2, xpc2, ypc2, euler, options);
end

