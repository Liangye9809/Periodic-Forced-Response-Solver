%**************************************************************************************************
% This file sloves the following equation using continuation method with Harmonic Balance Method:
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
%**************************************************************************************************
% INPUTS(getting from matlab Workspace):
% * CB.CBmods: structure contains the full-stuck elastic modes Phi and constrained modes Psi
% * CB.CB_MK: structure contains criag-bampton M K matrices
% * CB.CB_F: structure contains criag-bampton Fa Fx column vectors
% * xe0: preload displacement of elastic part
% * Rx: preload reaction forces act in contact part 
% 
%   OUTPUTS(save to .mat variables in current folder):
% * x_cont: solution of every ds in dof order with H harmonics fouriers' coefficient
% * omega_cont: corresponding omega of every ds
% 
% Written by Liu Liangye on June 03, 2025
% *************************************************************************************************


params.func.CBmods = CB.CBmods;
params.func.CB_MK = CB.CB_MK;
params.func.CB_F = CB.CB_F;
HBMstruct = struct('H', H, 'N', N, 'Nx', Nx, 'Na', Na, 'xi', xi, 'H_F_ext' ,H_F_ext, 'CB_F', params.func.CB_F);
params.func.HBM = HBM(HBMstruct); % H, N, Nc, Na, xi, Idx, E, EH
params.func.static.xp0 = xp0;
params.Newton = struct('epsx', epsx, 'epsf', epsf, 'maxiter', maxiter);
if ~exist('w', 'var')
    w = [0;0];
end
Coulombstruct = struct('kn', kn, 'xn0', xn0, 'mu', mu, 'kt', kt, 'Nx', Nx, 'w', w);
params.func.fc = CoulombFrictionParas(Coulombstruct); 
% fcvariables = CoulombFrictionParas(Coulombstruct); 
% params.func.fc.kn = fcvariables.kn;
% params.func.fc.xn0 = fcvariables.xn0;
% params.func.fc.mu = fcvariables.mu;
% params.func.fc.kt = fcvariables.kt;
% params.func.fc.w = fcvariables.w;
[xp, gxp, w] = Preloadxp(Rx, params); % preload displacement and preload forces in contact part
params.func.fc.w = w; % update w
params.func.static.preload = struct('xe0', xe0, 'Rx', Rx, 'xp', xp, 'gxp', gxp);
if size(x0,1) == 1
    x0 = x0 * ones((Na + 3 * Nx)*(2 * H + 1), 1);
end

tx0 = zeros(size(x0));

% define tlambda0 direction by lambda_end - lambda0
if omega_end - omega_0 < 0
    tomega0 = -1;
else
    tomega0 = 1;
end
params.cont = struct('ds', ds, 'maxstep', maxstep, 'omega_0', omega_0, 'omega_end', omega_end, ...
                     'x0', x0, 'tx0', tx0, 'tomega0', tomega0);
[x_cont, omega_cont] = continuation(@HBMFUNC, @HBMJACOB, @HBMJOmega, params);
