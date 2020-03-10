function euler_emsoft = converteuler_oim2emsoft(euler_oim, symID, lpangles)
%CONVERTEULER_OIM2EMSOFT Convert euler angles from OIM to EMsoft convention
%%% Inputs:
% -euler_oim: [phi1 PHI phi2] (deg) in OIM convention
% -symID: OIM symmetry ID from ang file
% -lpangles: lattice parameters [alpha beta gamma] (deg)
%%% Outputs:
% -euler_emsoft: [phi1 PHI phi2] (deg) in EMsoft convention

%%% Only monoclinic b (symID=20) is implemented so far
% Requires eu2om.m and om2eu.m
% 10/21/19

if symID==20
    Reuler = eu2om(euler_oim*pi/180);  % rotation matrix
    angle = -(lpangles(2)-90)*pi/180;   % beta (rad)
    Rrot = axang2rotm([0 1 0 angle]);   % rotation matrixfrom OIM to EMsoft
    R = Rrot*Reuler;    % perform rotation
    euler_emsoft = om2eu(R)*180/pi;     % Euler angles (deg) in EMsoft convention
elseif symID==1 || symID==2 || symID==3 || symID==32 || symID==6 || symID==62
    warning('Euler angles not explicitly converted to EMsoft convention, as this symmetry is not yet implemented. Results may be incorrect.');
    euler_emsoft = euler_oim;
else
    euler_emsoft = euler_oim;   % orthogonal crystal system, don't need to do anything
end


end

