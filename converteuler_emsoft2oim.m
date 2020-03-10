function euler_oim = converteuler_emsoft2oim(euler_emsoft, symID, lpangles)
%CONVERTEULER_EMSOFT2OIM Convert euler angles from EMsoft to OIM convention
%%% Inputs:
% -euler_emsoft: [phi1 PHI phi2] (deg) in OIM convention
% -symID: OIM symmetry ID from ang file
% -lpangles: lattice parameters [alpha beta gamma] (deg)
%%% Outputs:
% -euler_oim: [phi1 PHI phi2] (deg) in EMsoft convention

%%% Only monoclinic b (symID=20) is implemented so far
% Requires eu2om.m and om2eu.m
% 12/11/19

if symID==20
    Reuler = eu2om(euler_emsoft*pi/180);  % rotation matrix
    angle = (lpangles(2)-90)*pi/180;   % beta (rad)
    Rrot = axang2rotm([0 1 0 angle]);   % rotation matrix from EMsoft to OIM
    R = Rrot*Reuler;    % perform rotation
    euler_oim = om2eu(R)*180/pi;     % Euler angles (deg) in EMsoft convention
elseif symID==1 || symID==2 || symID==3 || symID==32 || symID==6 || symID==62
    warning('Euler angles not explicitly converted to EMsoft convention, as this symmetry is not yet implemented. Results may be incorrect.');
    euler_oim = euler_emsoft;
else
    euler_oim = euler_emsoft;   % orthogonal crystal system, don't need to do anything
end


end

