clear
clc
%% load data
Data
ReadFromCSV
Rx = 10*Rx;
% Nondimentionalization
%% preload
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
params.func.fc = CoulombFrictionParas(Coulombstruct); % handle classdef

[xp, gxp, w] = Preloadxp(Rx, params); % preload displacement and preload forces in contact part
params.func.fc.w = w;
params.func.static.preload = struct('xe0', xe0, 'Rx', Rx, 'xp', xp, 'gxp', gxp);
%%
M = [eye(5), CB.CB_MK.Max;
     CB.CB_MK.Max', CB.CB_MK.Mxx];
K = [diag(CB.CB_MK.Kaa), zeros(5,12);
     zeros(12,5), CB.CB_MK.Kxx];

F = [CB.CB_F.Fa;
     CB.CB_F.Fx];
% invM = inv(M);
F = M \ F;
F = [zeros(17,1); F];

% A_ = [zeros(17), eye(17);
%      -(invM * K), -xi * (invM * K)];

A = [zeros(17), eye(17);
     M \ -K, M \ (xi * -K)];

% R = [zeros(5,1); Rx];
% R = M \ R;
% R = [zeros(17,1); R];

y0 = zeros(17*2, 1);
y(:,1,1) = y0;
t(1,1) = 0;
i = 1;
omega_cont = [];
for omega = omega_0:omega_end
    % omega = omega * sqrt(CB.CB_MK.Kaa(1));
    omega_cont = [omega_cont, omega];
    T = 2*pi / omega;
    nstep = 300;
    nrelax = 5;
    dt = T / nstep;
    expdtAhalf = expm(0.5*dt*A);
    for k = 1:nrelax*nstep
        x = y(6:17, k, i)';
        Gstruct = g(x + xp', params.func);
        G = Gstruct.F - gxp';
        params.func.fc.w = Gstruct.w;
        G = [zeros(5,1); G'];
        G = M \ G;
        G = [zeros(17,1); G];
        yhalf = expdtAhalf * (y(:, k, i) + 0.5 * dt * (F * sin(omega * t(i, k)) - G));

        x = yhalf(6:17)';
        Gstruct = g(x + xp', params.func);
        G = Gstruct.F - gxp';
        params.func.fc.w = Gstruct.w;
        G = [zeros(5,1); G'];
        G = M \ G;
        G = [zeros(17,1); G];
        y(:, k+1, i) = expdtAhalf * (expdtAhalf * y(:, k, i) + dt * (F * sin(omega * (0.5*dt + t(i, k))) - G));

        t(i, k+1) = t(i, k) + dt;
    end
    if omega >= omega_end
        break;
    end
    i = i + 1;
    t(i, 1) = 0;
    y(:, 1, i) = y(:, end, i-1);
end

function G = phity(x, p)
        xp = p.static.preload.xp;
        gxp = p.static.preload.gxp;
        Gstruct = g(x + xp', p);
        G = Gstruct.F - gxp';
        params.func.fc.w = Gstruct.w;
        G = [zeros(5,1); G'];
        G = M \ G;
        G = [zeros(17,1); G];
end
%%
clear Amax
for j = 1:size(y,3) % omega number
    Amax(:, j) = max(abs(y(:,:,j)), [], 2);
end
omega_cont = omega_0:0.002:omega_end;
plot(omega_cont, Amax(1,:));