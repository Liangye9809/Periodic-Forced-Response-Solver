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

% omega02 = CB.CBmods.Kaa(1);
% alpha = abs(CB.CB_F.Fa(1)) / omega02;
% beta = abs(mu(1) * max(abs(Rx)) / kt(1));
% 
% 
% CB.CB_F.Fa = CB.CB_F.Fa / (alpha * omega02);
% CB.CB_F.Fx = CB.CB_F.Fx * beta / (alpha^2 * omega02);
% 
% Rx = Rx * beta / (alpha^2 * omega02);
% 
% CB.CB_MK.Max = CB.CB_MK.Max * (beta / alpha);
% CB.CB_MK.Mxx = CB.CB_MK.Mxx * (beta / alpha)^2;
% 
% CB.CB_MK.Kaa = CB.CB_MK.Kaa / omega02;
% CB.CBmods.Kaa = CB.CB_MK.Kaa;
% CB.CB_MK.Kxx = CB.CB_MK.Kxx * (beta / alpha)^2 / omega02;
% 
% kn = kn * (beta / alpha)^2 / omega02;
% kt = kt * (beta / alpha)^2 / omega02;
% 
% xi = xi * sqrt(omega02);
% 
% params.func.Nondimention.beta = beta;
% params.func.Nondimention.alpha = alpha;
% params.func.Nondimention.omega02 = omega02;
