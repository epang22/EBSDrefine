% Copyright (c) 2015, Marc De Graef
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


function q = om2eu(om)
% OM2EU Convert from rotation matrix to Euler angles
% From Marc De Graef 3Drotations github
% INPUT: 3x3 rotation matrix
% OUTPUT: [phi1 PHI phi2] in radians

thr = 1e-8;

if abs(abs(om(3,3))-1.0)>thr
    
    zeta = 1/sqrt(1-om(3,3)^2);
    q = [atan2(om(3,1)*zeta , -om(3,2)*zeta), acos(om(3,3)), ...
         atan2(om(1,3)*zeta , om(2,3)*zeta)];
    
else
  % we arbitrarily assign the entire angle to phi_1   
    if abs(om(3,3)-1.0)<thr
    
        q = [atan2(om(1,2),om(1,1)),0.0,0.0];
        
    else
    
        q = [-atan2(-om(1,2),om(1,1)),pi,0.0];
    end

end

% reduce Euler angles to definition ranges (and positive values only)
if (q(1)<0.0) 
    q(1) = mod(q(1)+100.0*pi,2.0*pi);
end

if (q(2)<0.0) 
    q(2) = mod(q(2)+100.0*pi,pi);
end

if (q(3)<0.0) 
    q(3) = mod(q(3)+100.0*pi,2.0*pi);
end