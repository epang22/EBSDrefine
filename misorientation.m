function misorient = misorientation(euler1, euler2)
%MISORIENTATION Compute misorientation between euler angles
% 2/7/20 (Edward Pang, MIT)
%
%%% Inputs:
% -euler1: [phi1 PHI phi2] (deg)
% -euler2: [phi1 PHI phi2] (deg)
%%% Outputs:
% -misorient: misorientation (deg)



% numerical tolerance threshold
thr = 1e-6;     % to deal with single precision error


% compute misorientation
q1 = eu2qu(euler1*pi/180);
q2 = eu2qu(euler2*pi/180);
qdot = dot(q1,q2);

% if qdot>1 (euler1 and euler2 are the same angle, numerical error)
if qdot>1
    if qdot-1<thr
        misorient = 0;
        return
    else
        % issue beyond numerical error
        warning('dp > 1+eps');
        % continue calcs, which will give complex numbers
    end
end

% if qdot~=1, continue...
misorient = 2*acosd(qdot);
misorient = min(misorient,360-misorient);   % put in range 0..180


