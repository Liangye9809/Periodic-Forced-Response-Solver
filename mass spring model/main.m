%% clear workspace and close all

clear
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


    
% dataname = 'Mesh32x32_CP' + string(Nx) + '.mat';
% % dataname = 'Mesh10x13_CP' + string(Nx) + '.mat';
% 
% tep = pwd;
% cd FEM/
% load(dataname);
% cd(tep)

Data




%% Craig-Bampton reduction from FEM matrices


CriagBamptonReduction


% load('CB.mat');
% load('Rx.mat');
% load('xe0.mat');



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
% omega_plot = 4100;

omega_0 = 0.1;
omega_end = 1.5;

% omega_0 = 0.2;
% omega_end = 2.3;

ContinuationCalculation


% t = [0:N-1]' * (2 * pi / N);
% xn = (- 4 * sin(sin(t)) + 1) .* 0.05; % separation to stick
% xt = (2 * exp(cos(t + 1)) - 3) .* 0.05; % separation to stick
% Xn = params.func.HBM.EH * xn;
% Xt = params.func.HBM.EH * xt;
% % figure;
% % plot(t, xt, 'b-'), hold on;
% % plot(t, xn, 'r-'), grid on;
% % legend('xt', 'xn');
% % for Na = 1
% params.func.HBM.fftfa = Xt / (alpha * omega02);
% Xt2 = zeros(2 * H + 1, 1);
% % Xt2(2) = 1 * 0.05;
% Xt2 = Xt;
% params.func.HBM.fftfx = 0.5 .* [Xt; Xt2; Xn] .* beta / (alpha^2 * omega02);

% Analytical
% [x_cont, omega_cont, k_cont, w_cont, stick_cont, slipP_cont, slipM_cont, gap_cont] = continuation(@HBMFUNC, @HBMJACOB, @HBMJOmega, params);

% Numerical
% [x_cont, omega_cont, k_cont] = continuation(@HBMFUNC, @HBMJACOB, @HBMJOmega, params);

% fixed Numerical
% [x_cont, omega_cont, k_cont, ~, ~] = continuation(@HBMFUNC, @HBMJACOB, @HBMJOmega, params);

CaseTime = toc;

%%

% HBMPostProcessing;
% HBMPostProcessing_Numerical;
% HBMPostProcessing_Numerical_fixed;


