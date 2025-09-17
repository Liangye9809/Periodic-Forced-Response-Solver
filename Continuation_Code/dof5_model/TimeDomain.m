clear
clc
%% load data
Data
CriagBamptonReduction

Nondimentionalization
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
M = [eye(2), CB.CB_MK.Max;
     CB.CB_MK.Max', CB.CB_MK.Mxx];
K = [diag(CB.CB_MK.Kaa), zeros(2,3);
     zeros(3,2), CB.CB_MK.Kxx];

F = [CB.CB_F.Fa;
     CB.CB_F.Fx];
% invM = inv(M);
F = M \ F;
F = [zeros(5,1); F];

% A_ = [zeros(17), eye(17);
%      -(invM * K), -xi * (invM * K)];

A = [zeros(5), eye(5);
     M \ (-K), M \ (xi * -K)];

% R = [zeros(5,1); Rx];
% R = M \ R;
% R = [zeros(17,1); R];

nstep = 200;
nrelax = 200;
nomega = 100;
y = zeros(5*2, 1+nstep*nrelax, nomega+1);
y0 = zeros(5*2, 1);
y(:,1,1) = y0;
t(1,1) = 0;
i = 1;
omega_cont = [];
domega = (omega_end - omega_0) / nomega;


%% Exponential Method with cubic

% tic;
% for omega = omega_0:domega:omega_end
%     omega_cont = [omega_cont, omega];
%     T = 2*pi / omega;
%     % nstep = 200;
%     % nrelax = 600;
%     dt = T / nstep;
%     R = expm(dt*A); 
%     for k = 1:nrelax*nstep
% 
%         x = y(3:5, k, i)';
%         Gstruct = g(x + xp', params.func);
%         G = Gstruct.F - gxp';
%         params.func.fc.w = Gstruct.w;
%         G = [zeros(2,1); G'];
%         G = M \ G;
%         G = [zeros(5,1); G];
% 
%         phi_ty = F * (sin(omega * t(i, k)) + sin(2 * omega * t(i, k))) - G;
%         y(:, k+1, i) = R * (y(:, k, i) + dt * phi_ty);
%         t(i, k+1) = t(i, k) + dt;
%     end
%     disp('omega ' + string(omega) + ' a1(end) ' + string(y(1, end, i)));
%     if omega >= omega_end
%         break;
%     end
%     i = i + 1;
%     t(i, 1) = 0;
%     y(:, 1, i) = y(:, end, i-1);
% 
% end
% 
% toc;


%% Exponential Method with friction
tic;
for omega = omega_0:domega:omega_end
    omega_cont = [omega_cont, omega];
    T = 2*pi / omega;
    % nstep = 200;
    % nrelax = 600;
    dt = T / nstep;
    R = expm(dt*A); 
    for k = 1:nrelax*nstep

        x = y(3:5, k, i)';
        Gstruct = gf(x + xp', params.func);
        G = Gstruct.F - gxp;
        params.func.fc.w = Gstruct.w;
        G = [zeros(2,1); G];
        G = M \ G;
        G = [zeros(5,1); G];

        % phi_ty = F * (sin(omega * t(i, k)) + sin(2 * omega * t(i, k))) - G;
        phi_ty = F * (sin(omega * t(i, k))) - G;
        y(:, k+1, i) = R * (y(:, k, i) + dt * phi_ty);
        t(i, k+1) = t(i, k) + dt;
    end
    disp('omega ' + string(omega) + ' a1(end) ' + string(y(1, end, i)));
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
plot(omega_cont, Amax(5,:),'ko');
grid on;
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
grid on;