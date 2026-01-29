%% clear workspace and close all
for H = 11:2:11
    for mainj = 8:8
clearvars -except mainj H
N = 2^mainj;
Nx = 1;
% for Nx = 4:4:36
% clearvars -except Nx

% clear
% close all
% clc
%**************************************************************************************************
% This file using continuation method to solve the equation under Harmonic Balance Method assumption :
% 
% |Mee Mec| |\ddot{e} |        |Kee Kec| |\dot{e} |   |Kee Kec| |e |   |  0  |   |Fe|        |Pe|
% |Mce Mcc|*|\ddot{xc}| + {xi}*|Kce Kcc|*|\dot{xc}| + |Kce Kcc|*|xc| + |g(xc)| = |Fc|*f(t) + |Pc|
% 
% Written by Liu Liangye on June 03, 2025
% *************************************************************************************************



%% Get parameters from file


    
dataname = 'Mesh32x32_CP' + string(Nx) + '.mat';
% dataname = 'Mesh10x13_CP' + string(Nx) + '.mat';

tep = pwd;
cd FEM/
load(dataname);
cd(tep)

Data




%% Craig-Bampton reduction from FEM matrices


CriagBamptonReduction

% load('CB.mat');
% load('Rx.mat');
% load('xe0.mat');

% ReadFromCSV % here we directly read the Criag-Bampton matrices from CSV file 

%% Dimensionless

%**************************************************************************************************
% This section make the dimensionless to the equation:
% 
% | I  Max| |\ddot{a} |        |Kaa  0 | |\dot{a} |   |Kaa  0 | |a |   |  0  |   |Fa|        | 0 |
% |Mxa Mxx|*|\ddot{xc}| + {xi}*| 0  Kxx|*|\dot{xc}| + | 0  Kxx|*|xc| + |g(xc)| = |Fx|*f(t) + |-Rx|
% 
% using the following changes:
% 
% t -> (1 / omega0) * t
% a -> alpha * a
% x -> beta * x
% 
% where:
% omega0 is selected to match that of the resonance of interest,
% alpha = fa / (omega0)^2
% beta = mu * max(abs(Rx)) / kt
% 
% Nondimentionalized Equation:
% |             I     (beta/alpha)Max| |\ddot{a} |               |(1/omega0)^2*Kaa                       0  | |\dot{a} |
% |(beta/alpha)Mxa (beta/alpha)^2*Mxx|*|\ddot{xc}| + {xi*omega0}*|              0  (beta/alpha/omega0)^2*Kxx|*|\dot{xc}|
% 
%    |(1/omega0)^2*Kaa                        0 | |a |   |               0           |   | (1/alpha/omega0^2)*Fa  |      
%  + |              0  (beta/alpha/omega0)^2*Kxx|*|xc| + |beta/(alpha*omega0)^2*g(xc)| = |(beta/alpha/omega0)^2*Fx|*f(t) 
% 
%                         | 0 |
% + beta/(alpha*omega0)^2*|-Rx|
%**************************************************************************************************
% INPUTS(getting from matlab Workspace):
% * CB.CBmods: structure contains the full-stuck elastic modes Phi and constrained modes Psi
% * CB.CB_MK: structure contains criag-bampton M K matrices
% * CB.CB_F: structure contains criag-bampton Fa Fx column vectors
% * Rx: preload reaction forces act in contact part
% 
%   OUTPUTS(save to .mat variables in current folder):
% * CB.CBmods: structure contains the full-stuck elastic modes Phi and constrained modes Psi
% * CB.CB_MK: structure contains criag-bampton M K matrices
% * CB.CB_F: structure contains criag-bampton Fa Fx column vectors
% * Rx: preload reaction forces act in contact part
% 
% Written by Liu Liangye on Sep 16, 2025
% *************************************************************************************************

Nondimensionalization
%% continuation calculation with HBM

tic;
omega_plot = 4100;
% switch Nx
%     case 4
%         omega_plot = 4210;
%     case 8
%         omega_plot = 4214.8;
%     case 12
%         omega_plot = 4216.2;
%     case 16
%         omega_plot = 4185;
%     case 20 
%         omega_plot = 4167.52;
%     case 24
%         omega_plot = 4175.41;
%     case 28
%         omega_plot = 4163.58;
%     case 32
%         omega_plot = 4170;
%     case 36
%         omega_plot = 4170;
% end

omega_0 = 3850 / sqrt(omega02);
% omega_end = 4400 / sqrt(omega02);
omega_end = omega_plot / sqrt(omega02);


ContinuationCalculation
CaseTime = toc;
% CaseInfo = 'Mesh32x32_Pe100each_Adof_CP' + string(Nx) + '_PreloadFixed_9points_H' + string(H) + '_N' + string(N);
% disp([CaseInfo, CaseTime]);


%%

HBMPostProcessing;


end
end

