function CB = cb_reduction_MK(FEMinput,CBmodes)
% *************************************************************************************************
% this function make the Braig-Bampton reduction of M and K matrices in the equations:
% 
% |Mee Mec| |\ddot{xe}|         |Kee Kec| |\dot{xe}|   |Kee Kec| |xe|   |  0  |   |F|
% |Mce Mcc|*|\ddot{xc}| + \{xi}*|Kce Kcc|*|\dot{xc}| + |Kce Kcc|*|xc| + |g(xc)| = |0|*f(t)
% 
% To:
% 
% | I  Max| |\ddot{a}|         |Kaa  0 | |\dot{a}|   |Kaa  0 | |a|   | 0  |   |Phi'*F|
% |Mxa Mxx|*|\ddot{x}| + \{xi}*| 0  Kxx|*|\dot{x}| + | 0  Kxx|*|x| + |g(x)| = |Psi'*F|*f(t)
% 
% where:
% 
% Max = Mxa' = Phi' * Mee * Psi + Phi' * Mec; 
% Mxx = Psi' * Mee * Psi + Psi' * Mec + Mce * Psi + Mcc; 
% Kxx = Kcc - Psi' * Kee * Psi; 
% Kaa = diag(omega^2);
% *************************************************************************************************
%   INPUTS:
% * FEMinput.Mee: mass matrix of elastic part in sparse form
% * FEMinput.Mec: mass matrix of coupling part in sparse form
% * FEMinput.Mcc: mass matrix of contact part in sparse form
% * FEMinput.Kee: stiffness matrix of elastic part in sparse form
% * FEMinput.Kcc: stiffness matrix of contact part in sparse form
% * CBmodes.Phi: first Na full-stuck elastic modes, Nc*Na matrix
% * CBmodes.Kaa: first Na frequencies square of full-stuck elastic modes, Na*1 vector
% * CBmodes.Psi: constrained modes, Ne*Nc matrix
%
%   OUTPUTS:
% * CB.Max: Na*Nc matrix in sparse form
% * CB.Mxa: inverse of Max
% * CB.Mxx: Nc*Nc matrix in sparse form
% * CB.Kxx: Nc*Nc matrix in sparse form
% * CB.Kaa: first Na frequencies square of full-stuck elastic modes, Na*1 vector
% 
% Written by Liu Liangye on April 07, 2025
% *************************************************************************************************

    Mee = FEMinput.Mee;
    Mec = FEMinput.Mec;
    Mcc = FEMinput.Mcc;
    Mce = Mec';
    Kee = FEMinput.Kee;
    Kcc = FEMinput.Kcc;

    Phi = CBmodes.Phi;
    Psi = CBmodes.Psi;
    if issymmetric(Mcc) && issymmetric(Kcc)
        Max = Phi'*Mee*Psi + Phi'*Mec; 
        Mxx = Psi'*Mee*Psi + Psi'*Mec + Mce*Psi + Mcc; 
        Mxa = Max';
        Kxx = Kcc - Psi'*Kee*Psi; 
        Kaa = CBmodes.Kaa;
        
        CB.Max = Max;
        CB.Mxa = Mxa;
        CB.Mxx = Mxx;
        CB.Kxx = Kxx;
        CB.Kaa = Kaa;
    else
        error('Matrix is not symmetric, please check');
    end
end


