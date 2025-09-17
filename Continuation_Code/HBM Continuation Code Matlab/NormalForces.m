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
 



function FN = NormalForces(xn, kn, xn0)
    u = xn - xn0;
    FN = max(0, kn .* u);
end