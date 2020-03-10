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


function q = eu2om(eu)
% EU2OM Convert from Euler angles to rotation matrix
% From Marc De Graef 3Drotations github
% INPUT: [phi1 PHI phi2] in radians
% OUTPUT: 3x3 rotation matrix

thr = 1e-10;

c1 = cos(eu(1));
c3 = cos(eu(3));
c2  = cos(eu(2));

s1 = sin(eu(1));
s3 = sin(eu(3));
s2  = sin(eu(2));

q = [ c1*c3-s1*c2*s3,  s1*c3+c1*c2*s3, s2*s3; ...
     -c1*s3-s1*c2*c3, -s1*s3+c1*c2*c3, s2*c3; ...
           s1*s2    ,      -c1*s2    ,  c2   ];

for i=1:3
  for j=1:3
    if (abs(q(i,j))< thr) 
        q(i,j) = 0.0;
    end
  end
end