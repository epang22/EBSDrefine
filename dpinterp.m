function [maxinterppoint, maxgridpoint, maxgridpointindices, F, exitflag] = ...
    dpinterp(dp_3d, test_angles_1Drodrigues, q1)
%DPINTERP
% Take dp vs orientation grid data and interpolate
%%% Inputs:
% -dp_3d: NxNxN array of dot products
% -test_angles_1Drodrigues: 1xm array of rodrigues grid in each dim about 0
% -q1: quaternion of grid center orientation
%%% Outputs:
% -maxinterppoint: 1x4 array containing phi1, PHI, phi2 (all in deg), dp of
%   interpolated max
% -maxgridppoint: 1x4 array containing phi1, PHI, phi2 (all in deg), dp of
%   best grid point
% -maxgridpointindices: 1x3 array containing 3D indices of max grid point
% -F: griddedInterpolant object
% -exitflag: status of interpolation. 0=ok, 1=failed (out of range 0..1), 
%   2=failed (worse than grid point), 3=warning out of domain range
% 12/27/19 (Edward Pang, MIT)



%%% Figure out max actual data point
[~,I] = max(reshape(dp_3d,1,[]));
[I1,I2,I3] = ind2sub(size(dp_3d),I);

% Compute euler angles of max grid point
x = test_angles_1Drodrigues(I1);    % Rodrigues vector components (grid centered about origin)
y = test_angles_1Drodrigues(I2);
z = test_angles_1Drodrigues(I3);
w = 1/sqrt(1+x^2+y^2+z^2);          % Calculate first quaternion entry
q2 = [w w*x w*y w*z];   % Quaternion of grid point (centered around north pole)
qrot = [q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)-q1(4)*q2(4), ...
    q1(1)*q2(2)+q1(2)*q2(1)+q1(3)*q2(4)-q1(4)*q2(3), ...
    q1(1)*q2(3)+q1(3)*q2(1)-q1(2)*q2(4)+q1(4)*q2(2), ...
    q1(1)*q2(4)+q1(4)*q2(1)+q1(2)*q2(3)-q1(3)*q2(2)];   % quaternion grid rotated to be centered on starting Euler angles
euler_maxgrid = qu2eu(qrot)*180/pi;     % Euler angles of grid point (in degrees)

% max grid point
maxgridpoint = [euler_maxgrid(1) euler_maxgrid(2) euler_maxgrid(3) dp_3d(I1,I2,I3)];
maxgridpointindices = [I1 I2 I3];


%%% Interpolate
F = griddedInterpolant({test_angles_1Drodrigues,test_angles_1Drodrigues,test_angles_1Drodrigues},dp_3d,'cubic');  % interpolate

% Define some quantities for optimization of interpolated surface
options = optimset('Display','off','TolFun',1e-8);  % Define optimization parameters (fminsearch)

% Find max interpolated point
fun = @(x)-F(x);    % Define function to minimize
x0 = [x y z];     % starting point of optimization
[testminparam,~] = fminsearch(fun,x0,options);

% Compute euler angles of max interpolated point
x = testminparam(1);    % Rodrigues vector components (grid centered about origin)
y = testminparam(2);
z = testminparam(3);
w = 1/sqrt(1+x^2+y^2+z^2);          % Calculate first quaternion entry
q2 = [w w*x w*y w*z];   % Quaternion of grid point (centered around north pole)
qrot = [q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)-q1(4)*q2(4), ...
    q1(1)*q2(2)+q1(2)*q2(1)+q1(3)*q2(4)-q1(4)*q2(3), ...
    q1(1)*q2(3)+q1(3)*q2(1)-q1(2)*q2(4)+q1(4)*q2(2), ...
    q1(1)*q2(4)+q1(4)*q2(1)+q1(2)*q2(3)-q1(3)*q2(2)];   % quaternion grid rotated to be centered on starting Euler angles
euler_maxinterp = qu2eu(qrot)*180/pi;     % Euler angles of max interp point (in degrees)

% max interpolated point
maxinterppoint = [euler_maxinterp(1) euler_maxinterp(2) euler_maxinterp(3) F([x y z])];


% check interpolation for errors
Rmax = max(test_angles_1Drodrigues);    % max value of Rodrigues vector component (+/-Rmax is domain bound)

if maxinterppoint(4) < 0 || maxinterppoint(4) > 1
    exitflag = 1;   % interpolation failed (out of range 0..1)
elseif maxinterppoint(4) < maxgridpoint(4)
    exitflag = 2;   % interpolation worse than grid point
elseif x > Rmax || x < -Rmax || y > Rmax || y < -Rmax || z > Rmax || z < -Rmax 
    exitflag = 3;   % interpolated max out of domain
else
    exitflag = 0;   % ok
end


