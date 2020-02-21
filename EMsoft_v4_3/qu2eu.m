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


function q = qu2eu(qq)
% QU2EU Convert from quaternions to Euler angles
% From Marc De Graef 3Drotations github
% INPUT: [q0 q1 q2 q3] unit quaternion
% OUTPUT: [phi1 PHI phi2] in radians

global epsijk

q0 = qq(1);
q1 = qq(2);
q2 = qq(3);
q3 = qq(4);

q03 = q0^2+q3^2;
q12 = q1^2+q2^2;
chi = sqrt(q03*q12);

if chi==0 && q12==0
    q = [atan2(-2*epsijk*q0*q3,q0^2-q3^2),0,0];
   
elseif chi==0 && q03==0
    q = [atan2(2*q1*q2,q1^2-q2^2),pi,0];
    
else 
    q = [atan2((q1*q3-epsijk*q0*q2)/chi,(-epsijk*q0*q1-q2*q3)/chi), ...
         atan2(2*chi,q03-q12),...
         atan2((epsijk*q0*q2+q1*q3)/chi,(-epsijk*q0*q1+q2*q3)/chi)];
end

% if (chi==0.0) 
%   if (q12==0.0) 
%    if (epsijk==1) 
%     Phi = 0.0;
%     phi2 = 0.0;                 % arbitrarily due to degeneracy
%     phi1 = atan2(-2.0*q0*q3,q2^2-q3^2);
%    else
%     Phi = 0.0;
%     phi2 = 0.0;                 % arbitrarily due to degeneracy
%     phi1 = atan2( 2.0*q0*q3,q2^2-q3^2);
%    end
%   else
%    if (epsijk==1) 
%     Phi = pi;
%     phi2 = 0.0;                 % arbitrarily due to degeneracy
%     phi1 = atan2(2.0*q1*q2,q1^2-q2^2);
%    else
%     Phi = pi;
%     phi2 = 0.0;                 % arbitrarily due to degeneracy
%     phi1 = atan2(2.0*q1*q2,q1^2-q2^2);
%    end
%   end
% else            % this is not a special degenerate case
%   if (epsijk==1) 
%     Phi = atan2( 2.0*chi, q03-q12 );
%     chi = 1/chi;
%     phi1 = atan2( (-q0*q2+q1*q3)*chi, (-q0*q1-q2*q3)*chi );
%     phi2 = atan2( (q0*q2+q1*q3)*chi, (-q0*q1+q2*q3)*chi );
%   else
%     Phi = atan2( 2.0*chi, q03-q12 );
%     chi = 1/chi;
%     phi1 = atan2( (q0*q2+q1*q3)*chi, (q0*q1-q2*q3)*chi );
%     phi2 = atan2( (-q0*q2+q1*q3)*chi, (q0*q1+q2*q3)*chi );
%   end
% end
% 
% q = [ phi1, Phi, phi2 ];

% reduce Euler angles to definition ranges (and positive values only)
if (q(1)<0.0) 
    q(1) = mod(q(1)+100.0*pi,2*pi);
end
if (q(2)<0.0) 
    q(2) = mod(q(2)+100.0*pi,pi);
end
if (q(3)<0.0) 
    q(3) = mod(q(3)+100.0*pi,2*pi);
end
