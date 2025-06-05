%**************************************************************************************************
% This file execute the Criag-Bampton reduction from equation:
% 
% |Mee Mec| |\ddot{e} |         |Kee Kec| |\dot{e} |   |Kee Kec| |e |   |  0  |   |Fe|        |Pe|
% |Mce Mcc|*|\ddot{xc}| + \{xi}*|Kce Kcc|*|\dot{xc}| + |Kce Kcc|*|xc| + |g(xc)| = |Fc|*f(t) + |Pc|
% 
% to equation:
% 
% | I  Max| |\ddot{a} |         |Kaa  0 | |\dot{a} |   |Kaa  0 | |a |   |  0  |   |Fa|        | 0 |
% |Mxa Mxx|*|\ddot{xc}| + \{xi}*| 0  Kxx|*|\dot{xc}| + | 0  Kxx|*|xc| + |g(xc)| = |Fx|*f(t) + |-Rx|
% 
% where:
%
% |Kee Kec| |xe0|   |Pe|   |0 |
% |Kce Kcc|*| 0 | = |Pc| + |Rx|
%
%  xe = e - xe0;
% 
% |xe|   |Phi Psi| |a |         |a |
% |xc| = | 0   I |*|xc| = Tcb * |xc|
% 
% Phi is full-stuck elastic modes
% Psi is constrained modes
%**************************************************************************************************
% INPUTS(getting from matlab Workspace):
% * Mee: mass matrix in elastic-elastic part
% * Mec: mass matrix in elastic-contact part
% * Mcc: mass matrix in contact-contact part
% * Kee: stiffness matrix in elastic-elastic part
% * Kec: stiffness matrix in elastic-contact part
% * Kcc: stiffness matrix in contact-contact part
% * Pe: preload forces in elastic part
% * Pc: preload forces in contact part
% * Fe: Amplitude of external forces in elastic part
% * Fc: Amplitude of external forces in contact part
% * Na: number of CB modes
% 
%   OUTPUTS(save to .mat variables in current folder):
% * CBmods: structure contains the full-stuck elastic modes Phi and constrained modes Psi
% * CB_MK: structure contains criag-bampton M K matrices
% * CB_F: structure contains criag-bampton Fa Fx column vectors
% * xe0: preload displacement of elastic part
% * Rx: preload reaction forces act in contact part
% 
% Written by Liu Liangye on June 03, 2025
% *************************************************************************************************


FEMinput = struct('Mee',Mee, 'Mec', Mec, 'Mcc', Mcc, 'Kee', Kee, 'Kec', Kec, 'Kcc', Kcc); 
F.Fe = Fe; 
F.Fc = Fc; 
CBmods = CBmodes(FEMinput, Na);
CB_MK = cb_reduction_MK(FEMinput, CBmods); 
CB_F = cb_reduction_vector(F, CBmods); 
xe0 = Kee \ Pe; 
Rx = Kec' * xe0 - Pc; 

save CBmods.mat CBmods
save CB_MK.mat CB_MK
save CB_F.mat CB_F
save xe0.mat xe0
save Rx.mat Rx