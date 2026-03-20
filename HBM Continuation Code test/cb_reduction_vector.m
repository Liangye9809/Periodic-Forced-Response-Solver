function CB = cb_reduction_vector(F,CBmodes)
% *************************************************************************************************
% this function make the Braig-Bampton reduction of column vector F in the equations:
% 
% |Mee Mec| |\ddot{xe}|         |Kee Kec| |\dot{xe}|   |Kee Kec| |xe|   |  0  |   |Fe|        | 0 |
% |Mce Mcc|*|\ddot{xc}| + \{xi}*|Kce Kcc|*|\dot{xc}| + |Kce Kcc|*|xc| + |g(xc)| = |Fc|*f(t) + |-Rx|
% 
% To:
% 
% | I  Max| |\ddot{a} |         |Kaa  0 | |\dot{a} |   |Kaa  0 | |a |   |  0  |   |Fa|        | 0 |
% |Mxa Mxx|*|\ddot{xc}| + \{xi}*| 0  Kxx|*|\dot{xc}| + | 0  Kxx|*|xc| + |g(xc)| = |Fx|*f(t) + |-Rx|
% 
% where:
% 
% Fa = Phi' * Fe;
% Fx = Psi' * Fe + Fc;
% *************************************************************************************************
%   INPUTS:
% * F: external forces coefficients in elastic part, Ne*1 vector
% * CBmodes.Phi: first Na full-stuck elastic modes, Nc*Na matrix
% * CBmodes.Psi: constrained modes, Ne*Nc matrix
%
%   OUTPUTS:
% * CB.Fa: Na*1 vector
% * CB.Fx: Nc*1 vector
%
%   where Na is the number of CB modes; 
%   Nc is the dof of contact part, = 3 * Nx
% 
% Written by Liu Liangye on April 07, 2025
% *************************************************************************************************
    
    Phi = CBmodes.Phi;
    Psi = CBmodes.Psi;
    
    Fe = F.Fe;
    Fc = F.Fc;

    CB.Fa = Phi' * Fe;
    CB.Fx = Psi' * Fe + Fc;
end


