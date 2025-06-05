function CB = cb_reduction_vector(F,CBmodes)
% *************************************************************************************************
% this function make the Braig-Bampton reduction of column vector F in the equations:
% 
% |Mee Mec| |\ddot{xe}|         |Kee Kec| |\dot{xe}|   |Kee Kec| |xe|   |  0  |   |F|
% |Mce Mcc|*|\ddot{xc}| + \{xi}*|Kce Kcc|*|\dot{xc}| + |Kce Kcc|*|xc| + |g(xc)| = |0|*f(t)
% 
% To:
% 
% | I  Max| |\ddot{a}|         |Kaa  0 | |\dot{a}|   |Kaa  0 | |a|   | 0  |   |Fa|
% |Mxa Mxx|*|\ddot{x}| + \{xi}*| 0  Kxx|*|\dot{x}| + | 0  Kxx|*|x| + |g(x)| = |Fx|*f(t)
% 
% where:
% 
% Fa = Phi' * F;
% Fx = Psi' * F;
% *************************************************************************************************
%   INPUTS:
% * F: external forces coefficients in elastic part, Ne*1 vector
% * CBmodes.Phi: first Na full-stuck elastic modes, Nc*Na matrix
% * CBmodes.Psi: constrained modes, Ne*Nc matrix
%
%   OUTPUTS:
% * CB.Fa: Na*1 vector
% * CB.Fx: Nx*1 vector
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


