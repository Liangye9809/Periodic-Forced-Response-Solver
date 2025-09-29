clear
% clc
%% load data
Data
ReadFromCSV
% Rx = 10*Rx; % cubic only
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

nstep = 50000;
nrelax = 500;
nomega = 20;
% y = zeros(17*2, 1+nstep*nrelax, nomega+1);

t(1,1) = 0;
i = 1;
omega_cont = [];
domega = (omega_end - omega_0) / nomega;

% %% Exponential Method without gx
% tic;
% i = 1;
% for omega = omega_0:domega:omega_end
%     omega_cont = [omega_cont, omega];
%     T = 2*pi / omega;
%     % nstep = 200;
%     % nrelax = 600;
%     dt = T / nstep;
%     R = expm(dt*A); 
%     for k = 1:nrelax*nstep
%         y(:, k+1, i) = R * (y(:, k, i) + dt * (F * sin(omega * t(i, k))));
%         t(i, k+1) = t(i, k) + dt;
%     end
%     disp('omega ' + string(omega) + ' a1(end) ' + string(y(1, end, i)));
%     if omega >= omega_end
%         break;
%     end
%     i = i + 1;
%     t(i, 1) = 0;
%     y(:, 1, i) = y(:, end, i-1);
% end
% 
% toc;
% 
% %% Exponential Method with cubic (nondimensionalization)
% 
% % stability
% J = finite_diff_jac(@(x) g(x, params.func).F, xp);
% J = [zeros(5,17);
%       zeros(12,5), J];
% J = -M \ J;
% Jg = [zeros(17,17), zeros(17,17);
%       J, zeros(17,17)];
% T = 2*pi / omega_0;
% for i = 0:5
%     nstep = 100 * 10^i;
%     dt = T / nstep;
%     R = dt * expm(dt*A) * Jg;
%     eigR = eig(R);
%     disp([nstep, dt, max(abs(real(eigR)))]);
% end
% 
% tic;
% i = 1;
% for omega = omega_0:domega:omega_end
%     omega_cont = [omega_cont, omega];
%     T = 2*pi / omega;
%     nstep = 200;
%     nrelax = 600;
%     dt = T / nstep;
%     R = expm(dt*A); 
%     for k = 1:nrelax*nstep
% 
%         x = y(6:17, k, i)';
%         Gstruct = g(x + xp', params.func);
%         G = Gstruct.F - gxp';
%         params.func.fc.w = Gstruct.w;
%         G = [zeros(5,1); G'];
%         G = M \ G;
%         G = [zeros(17,1); G];
% 
%         phi_ty = F * sin(omega * t(i, k)) - G;
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

%% Exponential Method with friction (nondimensionalization)
% J = finite_diff_jac(@(x) gf(x, params.func).F, xp);
% J = [zeros(5,17);
%       zeros(12,5), J];
% J = -M \ J;
% Jg = [zeros(17,17), zeros(17,17);
%       J, zeros(17,17)];

% tic;
% for omega = omega_0:domega:omega_end
%     omega_cont = [omega_cont, omega];
%     T = 2*pi / omega;
%     dt = T / nstep;
%     R = expm(dt*A); 
%     if i == 1
%         nrelax = 500;
%         y0 = zeros(17*2,nrelax*nstep+1);
%         cd .\data
%         load('y0');
%         y0(:,1,1) = y0_;
%         cd ..\
%         for k = 1:nrelax*nstep
% 
%             x = y0(6:17, k)';
%             Gstruct = gf(x + xp', params.func);
%             G = Gstruct.F - gxp;
%             params.func.fc.w = Gstruct.w;
%             G = [zeros(5,1); G];
%             G = M \ G;
%             G = [zeros(17,1); G];
% 
%             phi_ty = F * sin(omega * t(i, k)) - G;
%             y0(:, k+1) = R * (y0(:, k) + dt * phi_ty);
%             t(i, k+1) = t(i, k) + dt;
%         end
%         nrelax = 30;
%         y(:,:,1) = y0(:,end-nrelax*nstep:end);
%         disp('omega ' + string(omega) + ' a1(end) ' + string(y(1, end, i)));
%         i = i + 1;
%         t(i, 1) = 0;
%         y(:, 1, i) = y(:, end, i-1);
%     else
%         for k = 1:nrelax*nstep
% 
%             x = y(6:17, k, i)';
%             Gstruct = gf(x + xp', params.func);
%             G = Gstruct.F - gxp;
%             params.func.fc.w = Gstruct.w;
%             G = [zeros(5,1); G];
%             G = M \ G;
%             G = [zeros(17,1); G];
% 
%             phi_ty = F * sin(omega * t(i, k)) - G;
%             y(:, k+1, i) = R * (y(:, k, i) + dt * phi_ty);
%             t(i, k+1) = t(i, k) + dt;
%         end
%         disp('omega ' + string(omega) + ' a1(end) ' + string(y(1, end, i)));
%         if omega >= omega_end
%             break;
%         end
%         i = i + 1;
%         t(i, 1) = 0;
%         y(:, 1, i) = y(:, end, i-1);
%     end
% 
% end
% 
% toc;
% plot(y0(1,:),'b'), hold on;

%% Exponential Method with friction (nondimensionalization, with initial condition)
tic;
clear y t
t = 0;
nrelax = 500;
y = zeros(17*2, nrelax*nstep+1);
cd .\data
load('y0');
y(:,1) = y0_;
cd ..\
%%
for omega = (omega_0+domega*2):domega:omega_end
    omega_cont = [omega_cont, omega];
    T = 2*pi / omega;
    dt = T / nstep;
    R = expm(dt*A); 

    

    for k = 1:nrelax*nstep

        x = y(6:17, k)';
        Gstruct = gf(x + xp', params.func);
        G = Gstruct.F - gxp;
        params.func.fc.w = Gstruct.w;
        G = [zeros(5,1); G];
        G = M \ G;
        G = [zeros(17,1); G];

        phi_ty = F * sin(omega * t) - G;
        y(:, k+1) = R * (y(:, k) + dt * phi_ty);
        t = t + dt;
    end
    disp('omega ' + string(omega) + ' a1(end) ' + string(y(1, end)));
    Amax(:, i) = max(abs(y(:,end-10*nstep:end)), [], 2);
    if omega >= omega_end
        break;
    end
    i = i + 1;
    t = 0;
    y(:, 1) = y(:, end);
    

end

toc;
% plot(y0(1,:),'b'), hold on;

% %% Implicit Euler Method without gx
% tic;
% for omega = omega_0:domega:omega_end
%     omega_cont = [omega_cont, omega];
%     T = 2*pi / omega;
%     % nstep = 200;
%     % nrelax = 500;
%     dt = T / nstep;
%     R = expm(dt*A); 
%     for k = 1:nrelax*nstep
%         y(:, k+1, i) = R * y(:, k, i) + dt * (F * sin(omega * t(i, k)));
%         t(i, k+1) = t(i, k) + dt;
%     end
%     if omega >= omega_end
%         break;
%     end
%     i = i + 1;
%     t(i, 1) = 0;
%     y(:, 1, i) = y(:, end, i-1);
% end
% 
% toc;
% 
% %% Implicit Euler Method with gx iteration
% tic;
% for omega = omega_0:domega:omega_end
%     omega_cont = [omega_cont, omega];
%     T = 2*pi / omega;
%     % nstep = 1000;
%     % nrelax = 100;
%     dt = T / nstep;
%     % expdtAhalf = expm(0.5*dt*A);
%     R = inv(eye(34) - dt*A);
%     for k = 1:nrelax*nstep
%         x = y(6:17, k, i)';
%         Gstruct = g(x + xp', params.func);
%         G = Gstruct.F - gxp';
%         params.func.fc.w = Gstruct.w;
%         G = [zeros(5,1); G'];
%         G = M \ G;
%         G = [zeros(17,1); G];
%         yp1 = R * (y(:, k, i) + dt * (F * sin(omega * (dt + t(i, k))) - G));
% 
%         % function iteration
%         % x = yp1(6:17)'; % phi(t,y+1)
%         % for j = 2:50
%         %     Gstruct = g(x + xp', params.func);
%         %     G = Gstruct.F - gxp';
%         %     params.func.fc.w = Gstruct.w;
%         %     G = [zeros(5,1); G'];
%         %     G = M \ G;
%         %     G = [zeros(17,1); G];
%         %     yp12 = R * (y(:, k, i) + dt * (F * sin(omega * (dt + t(i, k))) - G));
%         %     errory = norm(yp12 - yp1) / norm(yp1);
%         %     if max(abs(errory)) < 1e-3
%         %         break;
%         %     end
%         %     if j>49
%         %         warning('iteration cannot converge in omega %f of nstep %d',omega,k);
%         %     end
%         %     yp1 = yp12;
%         %     x = yp12(6:17)';
%         % end
% 
%         y(:, k+1, i) = yp1;
%         t(i, k+1) = t(i, k) + dt;
%     end
%     disp('omega ' + string(omega) + ' a1(end) ' + string(yp1(1)));
%     if omega >= omega_end
%         break;
%     end
%     i = i + 1;
%     t(i, 1) = 0;
%     y(:, 1, i) = y(:, end, i-1);
% end
% 
% toc;
% 
% %% RK2 method without gx
% tic;
% domega = 5;
% for omega = omega_0:domega:omega_end
%     omega_cont = [omega_cont, omega];
%     T = 2*pi / omega;
%     % nstep = 14500;
%     % nrelax = 50;
%     dt = T / nstep;
%     R = eye(34) + (dt/2)*A; % half step of Explicit Euler Method
%     for k = 1:nrelax*nstep
%         Ft = F * sin(omega * (dt + t(i, k)));
%         yhalf = R * y(:, k, i) + (dt/2) * Ft;
%         Ft = F * sin(omega * (dt/2 + t(i, k)));
%         y(:, k+1, i) = y(:, k, i) + dt * (A * yhalf + Ft);
%         t(i, k+1) = t(i, k) + dt;
%     end
%     if omega >= omega_end
%         break;
%     end
%     i = i + 1;
%     t(i, 1) = 0;
%     y(:, 1, i) = y(:, end, i-1);
% end
% 
% toc;

%%
% clear Amax
% for j = 1:size(y,3) % omega number
%     Amax(:, j) = max(abs(y(:,end-10*nstep:end,j)), [], 2);
% end
% omega_cont = omega_0:domega:omega_end;
figure
plot(omega_cont, Amax(1,:),'ko');
grid on;
%% oscilation over time
% figure;
% i = 1;
% tt = 0;
% for omega = omega_0:domega:omega_end
%     T = 2 * pi / omega;
%     dt = T / nstep;
%     tt = tt(end) + dt*[0:nstep*nrelax];
%     plot(tt, y(1,:,i)), hold on;
%     i = i + 1;
% end
% grid on;

%% oscilation over cycle number
% figure;
% i = 1;
% ss = 0;
% for omega = omega_0:domega:omega_end
%     T = 2 * pi / omega;
%     dt = T / nstep;
%     ss = ss(end) + [0:nstep*nrelax] ./ nstep;
%     plot(ss, y(1,:,i)), hold on;
%     i = i + 1;
% end
% grid on;
% xticks(0:30:630);
% xlabel('cycles');
% ylabel('a_1');
%%
figure;
ss = 0;


ss = ss(end) + [0:nstep*nrelax] ./ nstep;
plot(ss, y(1,:)), hold on;

grid on;
% xticks(0:30:630);
xlabel('cycles');
ylabel('a_1');