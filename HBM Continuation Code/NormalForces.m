% *************************************************************************************************
% this function calculate the normal forces of each contact points
% if normal displacement xn > lift-off displacement xn0
% FN = kn * (xn - xn0)
% else 
% FN = 0;
% *************************************************************************************************
%   INPUTS:
% * xn: normal displacement in N interval of time
% * kn: normal stiffness
% * xn0: lift-off displacement
% 
%   OUTPUTS:
% * FN: normal forces
% 
% Written by Liu Liangye on April 30, 2025
% *************************************************************************************************
 

% function FN = NormalForces(xn, kn, xn0)
%     u = xn - xn0;
%     if u < 0
%         FN = 0;
%     else
%         FN = kn .* u;
%     end
% end

% assume xn, kn and xn0 are all row vectors
% function FN = NormalForces(xn, kn, xn0)
%     u = xn - xn0;
%     u(u<0) = 0;
%     FN = kn .* u;
% end

function FN = NormalForces(xn, kn, xn0)
    u = xn - xn0;
    FN = max(0, kn .* u);
end