clear
clc
%% load data
Data
ReadFromCSV
Rx = 0*Rx;
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
params.func.fc = CoulombFrictionParas(Coulombstruct); 

[xp, gxp, w] = Preloadxp(Rx, params); % preload displacement and preload forces in contact part
params.func.fc.w = w; % update w
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
     M \ (-K), M \ (xi * -K)];

% R = [zeros(5,1); Rx];
% R = M \ R;
% R = [zeros(17,1); R];

y0 = zeros(17*2, 1);
y(:,1,1) = y0;
t(1,1) = 0;
i = 1;
omega_cont = [];
domega = (omega_end - omega_0) / 25;
%% 
tic;
for omega = omega_0:domega:omega_end
    % omega = omega * sqrt(CB.CB_MK.Kaa(1));
    omega_cont = [omega_cont, omega];
    T = 2*pi / omega;
    nstep = 300;
    nrelax = 100;
    dt = T / nstep;
    % expdtAhalf = expm(0.5*dt*A);
    % R = inv(eye(34) - dt*A);
    R = eye(34) - dt*A; % Implicit Euler Method
    for k = 1:nrelax*nstep
        % x = y(6:17, k, i)';
        % Gstruct = g(x + xp', params.func);
        % G = Gstruct.F - gxp';
        % params.func.fc.w = Gstruct.w;
        % G = [zeros(5,1); G'];
        % G = M \ G;
        % G = [zeros(17,1); G];
        % yp1 = R * y(:, k, i) + dt * (F * sin(omega * (dt + t(i, k))) - G);
        y(:, k+1, i) = R \ (y(:, k, i) + dt * (F * sin(omega * (dt + t(i, k)))));
        % x = yp1(6:17)'; % phi(t,y+1)
        % for j = 2:50
        %     Gstruct = g(x + xp', params.func);
        %     G = Gstruct.F - gxp';
        %     params.func.fc.w = Gstruct.w;
        %     G = [zeros(5,1); G'];
        %     G = M \ G;
        %     G = [zeros(17,1); G];
        %     yp12 = R * y(:, k, i) + dt * (F * sin(omega * (dt + t(i, k))) - G);
        %     errory = norm(yp12 - yp1) / norm(yp1);
        %     if max(abs(errory)) < 1e-3
        %         break;
        %     end
        %     if j>49
        %         warning('iteration cannot converge in omega %f of nstep %d',omega,k);
        %     end
        %     yp1 = yp12;
        %     x = yp12(6:17)';
        % end
        % y(:, k+1, i) = yp12;
        t(i, k+1) = t(i, k) + dt;
    end
    if omega >= omega_end
        break;
    end
    i = i + 1;
    t(i, 1) = 0;
    y(:, 1, i) = y(:, end, i-1);
end

toc;

%% RK2 method without gx
tic;
domega = 5;
for omega = omega_0:domega:omega_end
    omega_cont = [omega_cont, omega];
    T = 2*pi / omega;
    nstep = 14500;
    nrelax = 50;
    dt = T / nstep;
    R = eye(34) + (dt/2)*A; % half step of Explicit Euler Method
    for k = 1:nrelax*nstep
        Ft = F * sin(omega * (dt + t(i, k)));
        yhalf = R * y(:, k, i) + (dt/2) * Ft;
        Ft = F * sin(omega * (dt/2 + t(i, k)));
        y(:, k+1, i) = y(:, k, i) + dt * (A * yhalf + Ft);
        t(i, k+1) = t(i, k) + dt;
    end
    if omega >= omega_end
        break;
    end
    i = i + 1;
    t(i, 1) = 0;
    y(:, 1, i) = y(:, end, i-1);
end

toc;

%%
clear Amax
for j = 1:size(y,3) % omega number
    Amax(:, j) = max(abs(y(:,end-10*nstep:end,j)), [], 2);
end
% omega_cont = omega_0:domega:omega_end;
figure
plot(omega_cont, Amax(1,:),'ko');
%%
figure;
i = 1;
tt = 0;
for omega = omega_0:domega:omega_end
    T = 2 * pi / omega;
    dt = T / nstep;
    tt = tt(end) + dt*[0:nstep*nrelax];
    plot(tt, y(1,:,i)), hold on;
    i = i + 1;
end