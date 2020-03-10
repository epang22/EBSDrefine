function [test_angles, test_angles_1Drodrigues, quaternions, misorientations] = createorigrid(delta_angles, N_angles, euler_center)
%CREATEORIGRID Create test grid of angles
% First create grid about origin in Rodrigues space, then convert to
% quaternion space centered about north pole [1,0,0,0], then rotate to
% be centered around starting Euler angles and convert to Euler angles
% Inputs:
% -delta_angles: half-width of grid in deg
% -N_angles: number of points from euler to euler+delta_angles (grid will
% be 2N_angles+1 in each dim)
% -euler: [phi1 PHI phi2] in deg
% Outputs: 
% -test_angles: (2N+1)^3 by 3 array, each row is euler angles of a grid point in deg
% -test_angles_1Drodrigues: 1 by 2N+1 array, rodrigues vector components
% -quaternions: (2N+1)^3 by 4 array, each row is unit quaternion of a grid point
% -misorientations: (2N+1)^3 by 1 array, misorientation of each grid point (deg) from euler_center
% ]4/20/19 (Edward Pang, MIT)


if N_angles>0
    test_angles_1Drodrigues = -tand(delta_angles/2):tand(delta_angles/2)/N_angles:tand(delta_angles/2);
    test_angles = zeros((2*N_angles+1)^3,3);    % Initialize; each row is an angle
    quaternions = zeros((2*N_angles+1)^3,4);    % Initialize; each row is an angle
    misorientations = zeros((2*N_angles+1)^3,1);    % Initialize; columns are Euler angles;
    q1 = eu2qu(euler_center*pi/180);     % Quaternion of starting Euler angles (center of grid)
    for ii=1:2*N_angles+1
        for jj=1:2*N_angles+1
            for kk=1:2*N_angles+1
                x = test_angles_1Drodrigues(ii);    % Rodrigues vector components (grid centered about origin)
                y = test_angles_1Drodrigues(jj);
                z = test_angles_1Drodrigues(kk);
                w = 1/sqrt(1+x^2+y^2+z^2);          % Calculate first quaternion entry
                q2 = [w w*x w*y w*z];   % Quaternion of grid point (centered around north pole)
                qrot = [q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)-q1(4)*q2(4), ...
                    q1(1)*q2(2)+q1(2)*q2(1)+q1(3)*q2(4)-q1(4)*q2(3), ...
                    q1(1)*q2(3)+q1(3)*q2(1)-q1(2)*q2(4)+q1(4)*q2(2), ...
                    q1(1)*q2(4)+q1(4)*q2(1)+q1(2)*q2(3)-q1(3)*q2(2)];   % quaternion grid rotated to be centered on starting Euler angles
                euler = qu2eu(qrot)*180/pi;     % Euler angles of grid point (in degrees)
                
                test_angles((ii-1)*(2*N_angles+1)^2+(jj-1)*(2*N_angles+1)+kk,1) = euler(1);
                test_angles((ii-1)*(2*N_angles+1)^2+(jj-1)*(2*N_angles+1)+kk,2) = euler(2);
                test_angles((ii-1)*(2*N_angles+1)^2+(jj-1)*(2*N_angles+1)+kk,3) = euler(3);
                
                quaternions((ii-1)*(2*N_angles+1)^2+(jj-1)*(2*N_angles+1)+kk,:) = qrot;
                
                misorientations((ii-1)*(2*N_angles+1)^2+(jj-1)*(2*N_angles+1)+kk) = 2*acosd(dot(q1,qrot));
            end
        end
    end   
end

end

