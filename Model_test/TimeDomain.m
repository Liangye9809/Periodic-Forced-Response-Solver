clear
% clc
%% load data
Data

if ~exist('w', 'var')
    w = [0;0];
end
Coulombstruct = struct('kn', kn, 'xn0', xn0, 'mu', mu, 'kt', kt, 'Nx', Nx, 'w', w);
params.func.fc = CoulombFrictionParas(Coulombstruct); 

M = readmatrix('M.csv');
K = readmatrix('K.csv');
f = readmatrix('f.csv');
F = 0.3 * f;
gxp = readmatrix('gp.csv');
xp = readmatrix('xp.csv');

params.func.static.preload = struct('xp', xp, 'gxp', gxp);
%%

F = M \ F;
F = [zeros(51,1); F];

A = [zeros(51), eye(51);
     M \ (-K), M \ (xi * -K)];

nstep = 10000;
nrelax = 30;
nomega = 20;
y = zeros(51*2, 1+nstep*nrelax, nomega+1);

t(1,1) = 0;
i = 1;
omega_cont = [];
domega = (omega_end - omega_0) / nomega;

%% Exponential Method with cubic (nondimensionalization)

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

% stability
% J = finite_diff_jac(@(x) gf(x, params.func).F, xp');
% J = [zeros(3,51);
%       zeros(48,3), J];
% J = -M \ J;
% % Jg = [zeros(51,51), zeros(51,51);
% %       zeros(51,51), J];
% Jg = [zeros(51,51), zeros(51,51);
%       J, zeros(51,51)];
% T = 2*pi / omega_0;
% for i = 0:5
%     nstep = 100 * 10^i;
%     dt = T / nstep;
%     R = dt * expm(dt*A) * Jg;
%     eigR = eig(R);
%     disp([nstep, dt, max(abs(real(eigR)))]);
% end

tic;
for omega = omega_0:domega:omega_end
    omega_cont = [omega_cont, omega];
    T = 2*pi / omega;
    dt = T / nstep;
    R = expm(dt*A); 
    if i == 1
        nrelax = 500;
        y0 = zeros(51*2,nrelax*nstep+1);
        for k = 1:nrelax*nstep
    
            x = y0(4:51, k)';
            Gstruct = gf(x + xp', params.func);
            G = Gstruct.F - gxp;
            params.func.fc.w = Gstruct.w;
            G = [zeros(3,1); G];
            G = M \ G;
            G = [zeros(51,1); G];
    
            phi_ty = F * sin(omega * t(i, k)) - G;
            y0(:, k+1) = R * (y0(:, k) + dt * phi_ty);
            t(i, k+1) = t(i, k) + dt;
        end
        nrelax = 30;
        y(:,:,1) = y0(:,end-nrelax*nstep:end);
        disp('omega ' + string(omega) + ' a1(end) ' + string(y(1, end, i)));
        i = i + 1;
        t(i, 1) = 0;
        y(:, 1, i) = y(:, end, i-1);
    else
        for k = 1:nrelax*nstep
    
            x = y(4:51, k, i)';
            Gstruct = gf(x + xp', params.func);
            G = Gstruct.F - gxp;
            params.func.fc.w = Gstruct.w;
            G = [zeros(3,1); G];
            G = M \ G;
            G = [zeros(51,1); G];
    
            phi_ty = F * sin(omega * t(i, k)) - G;
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
grid on;
%% oscilation over time
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

%% oscilation over cycle number
figure;
i = 1;
ss = 0;
for omega = omega_0:domega:omega_end
    T = 2 * pi / omega;
    dt = T / nstep;
    ss = ss(end) + [0:nstep*nrelax] ./ nstep;
    plot(ss, y(1,:,i)), hold on;
    i = i + 1;
end
grid on;
xticks(0:30:630);
xlabel('cycles');
ylabel('a_1');