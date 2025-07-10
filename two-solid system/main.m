%% clear workspace and close all
clear
close all
clc
%**************************************************************************************************
% This file using continuation method with Harmonic Balance Method assumption to slove the equation:
% 
% |Mee Mec| |\ddot{e} |         |Kee Kec| |\dot{e} |   |Kee Kec| |e |   |  0  |   |Fe|        |Pe|
% |Mce Mcc|*|\ddot{xc}| + \{xi}*|Kce Kcc|*|\dot{xc}| + |Kce Kcc|*|xc| + |g(xc)| = |Fc|*f(t) + |Pc|
% 
% Written by Liu Liangye on June 03, 2025
% *************************************************************************************************



%% Get parameters from file

Data
% Biud_FEM_structure

%% Criag-Bampton reduction from FEM matrices

%**************************************************************************************************
% This section execute the Criag-Bampton reduction from equation:
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
% * FEM.Mee: mass matrix in elastic-elastic part
% * FEM.Mec: mass matrix in elastic-contact part
% * FEM.Mcc: mass matrix in contact-contact part
% * FEM.Kee: stiffness matrix in elastic-elastic part
% * FEM.Kec: stiffness matrix in elastic-contact part
% * FEM.Kcc: stiffness matrix in contact-contact part
% * FEM.Pe: preload forces in elastic part
% * FEM.Pc: preload forces in contact part
% * FEM.Fe: Amplitude of external forces in elastic part
% * FEM.Fc: Amplitude of external forces in contact part
% * FEM.Na: number of CB modes
% 
%   OUTPUTS(save to .mat variables in current folder):
% * CB.CBmods: structure contains the full-stuck elastic modes Phi and constrained modes Psi
% * CB.CB_MK: structure contains criag-bampton M K matrices
% * CB.CB_F: structure contains criag-bampton Fa Fx column vectors
% * xe0: preload displacement of elastic part
% * Rx: preload reaction forces act in contact part
% 
% Written by Liu Liangye on June 03, 2025
% *************************************************************************************************

% CriagBamptonReduction
ReadFromCSV
% Rx = 10*Rx; % cubic only
%% Dimensionless

%**************************************************************************************************
% This section make the dimensionless to the equation:
% 
% | I  Max| |\ddot{a} |         |Kaa  0 | |\dot{a} |   |Kaa  0 | |a |   |  0  |   |Fa|        | 0 |
% |Mxa Mxx|*|\ddot{xc}| + \{xi}*| 0  Kxx|*|\dot{xc}| + | 0  Kxx|*|xc| + |g(xc)| = |Fx|*f(t) + |-Rx|
% 
% to equation:
% 
% 
% 
% 
% where:
%
% 
% 
%
%  
% 
% 
% 
% 
% 
% 
%**************************************************************************************************
% INPUTS(getting from matlab Workspace):
% * 
% * 
% * 
% * 
% * 
% * 
% * 
% * 
% * 
% * 
% * 
% 
%   OUTPUTS(save to .mat variables in current folder):
% * 
% * 
% * 
% * 
% * 
% 
% Written by Liu Liangye on June 03, 2025
% *************************************************************************************************

Nondimentionalization
%% continuation calculation with HBM

%**************************************************************************************************
% This section using continuation method with Harmonic Balance Method assumption to slove the equation:
% 
% | I  Max| |\ddot{a}|         |Kaa  0 | |\dot{a}|   |Kaa  0 | |a|   |         0       |   |Fa|     
% |Mxa Mxx|*|\ddot{x}| + \{xi}*| 0  Kxx|*|\dot{x}| + | 0  Kxx|*|x| + |g(x + xp) - g(xp)| = |Fx|*f(t)
% 
% from:
% 
% | I  Max| |\ddot{a} |         |Kaa  0 | |\dot{a} |   |Kaa  0 | |a |   |  0  |   |Fa|        | 0 |
% |Mxa Mxx|*|\ddot{xc}| + \{xi}*| 0  Kxx|*|\dot{xc}| + | 0  Kxx|*|xc| + |g(xc)| = |Fx|*f(t) + |-Rx|
% 
% where:
%
% |Kaa  0 | |0 |   |  0  |   | 0 |
% | 0  Kxx|*|xp| + |g(xp)| = |-Rx|
%
% x = xc - xp
%**************************************************************************************************
% INPUTS(getting from matlab Workspace):
% * CB.CBmods: structure contains the full-stuck elastic modes Phi and constrained modes Psi
% * CB.CB_MK: structure contains criag-bampton M K matrices
% * CB.CB_F: structure contains criag-bampton Fa Fx column vectors
% * CB.xe0: preload displacement of elastic part
% * CB.Rx: preload reaction forces act in contact part 
% 
%   OUTPUTS(save to .mat variables in current folder):
% * x_cont: solution of every ds in dof order with H harmonics fouriers' coefficient
% * omega_cont: corresponding omega of every ds
% 
% Written by Liu Liangye on June 03, 2025
% *************************************************************************************************
tic;
profile on
ContinuationCaluculation
profile off
toc;
p = profile('info');
save('ProfileInfo.mat','p');
%%

HBMPostProcessing;

