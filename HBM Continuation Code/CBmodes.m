function CBmods = CBmodes(FEMinput, Na)
% *************************************************************************************************
% this function calculate the full-stuck elastic modes and constrain modes of the equations of motion:
% 
% |Mee Mec| |\ddot{e} |         |Kee Kec| |\dot{e} |   |Kee Kec| |e |   |  0  |   |Fe|        |Pe|
% |Mce Mcc|*|\ddot{xc}| + \{xi}*|Kce Kcc|*|\dot{xc}| + |Kce Kcc|*|xc| + |g(xc)| = |Fc|*f(t) + |Pc|
% 
% full-stuck elastic modes Phi are calculated by eigenvalue problem:
% (-omega^2 * Mee + Kee) * Phi = 0;
% 
% constrained modes Psi are calculated by static equilibrium equation:
% Kee * e + Kec * xc = 0;
% e = -Kee \ Kec * xc = Psi * xc;
% *************************************************************************************************
%   INPUTS:
% * FEMinput.Mee: mass matrix of elastic part in sparse form
% * FEMinput.Kee: stiffness matrix of elastic part in sparse form
% * FEMinput.Kec: stiffness matrix of coupling part in sparse form
% * Na: number of full-stuck elastic modes
% 
%   OUTPUTS:
% * CBmods.Phi: first Na full-stuck elastic modes, Nc*Na matrix
% * CBmods.Kaa: first Na frequencies square of full-stuck elastic modes, Na*1 vector
% * CBmods.Psi: constrained modes, Ne*Nc sparse matrix
% 
%   where Na is the number of CB modes; 
%   Nc is the dof of contact part, = 3 * Nx
%
% Written by Liu Liangye on April 07, 2025
% *************************************************************************************************
 
Mee = FEMinput.Mee;
Kee = FEMinput.Kee;
Kec = FEMinput.Kec;

if issymmetric(Mee) && issymmetric(Kee)
    [Phi, D] = eigs(Kee, Mee, Na, 'smallestabs');
    Kaa = diag(D);
    Psi = -Kee \ Kec;
    
    CBmods.Phi = Phi;
    CBmods.Psi = Psi;
    CBmods.Kaa = Kaa;
else
    error('Matrix is not symmetric, please check');
   
end


end

