clear
clc
%% for single dof of variable
eps = zeros(11, 10);
h = 10^(-8);
order = 1;
h_con = zeros(11, 10);
for iH = 1:11
    for iN = 1:12
        % N = 2 ^ (iN + 4);
        N = 256;
        H = iH;
        ih = iN + 0;
        h = 10^(-ih);
        dt = 2 * pi / N;
        t = 0:dt:(2 * pi - dt);
        xt = exp(sin(t))';
        dgt = dg(xt);
        % plot(t,x);
        % plot(t,dgt);
        
        % Numerical Jacobian
        [E, EH] = HBM.fft_matrices(N, H);
        X = EH * xt;
        JNL_N = finite_diff_jac(@(x) fftgx(x, E, EH), X, h, order);
        
        % Analytiacl Jacobian
        [E, EH] = HBM.fft_matrices(N, 2 * H);
        dGmn = EH * dgt;
        JNL_A = HBMJACOB_analytical(dGmn, H);

        eps(iH, iN) = norm(JNL_A - JNL_N) / norm(JNL_A);
        h_con(iH, iN) = h;
    end
end

%% for multi dof of variables
clear
clc
eps = zeros(11, 10);
h = 10^(-8);
order = 2;
h_con = zeros(11, 10);
N = 256;
dt = 2 * pi / N;
t = 0:dt:(2 * pi - dt);
xt = fxt(t);
dgt = dg(xt);
% plot(t,xt);
% plot(t,dgt);
for iH = 1:20
    for iN = 1:12
        
        H = iH;
        ih = iN + 0;
        h = 10^(-ih);

        % Numerical Jacobian
        [E, EH] = HBM.fft_matrices(N, H);
        X = EH * xt;
        X = X(:);
        % [F, xtF] = fftgx(X, E, EH); % for testing
        JNL_N = finite_diff_jac(@(x) fftgx(x, E, EH), X, h, order);
        
        % Analytiacl Jacobian
        [E, EH] = HBM.fft_matrices(N, 2 * H);
        dGmn = []; % initialize the size
        for iG = 1:3
            for jG = 1:3
                dGtemp = [];
                dGtemp(:, 1) = dgt(iG, jG, :);
                dGmn(iG, jG, :) = EH * dGtemp;
            end
        end
        JNL_A = HBMJACOB_analytical_muilti(dGmn, H);

        eps(iH, iN) = norm(JNL_A - JNL_N) / norm(JNL_A);
        h_con(iH, iN) = h;
    end
end

%% Coulomb friction
clear
clc
load('Pe100eachData_PeFixed_Omega_4160_Nx_1_H11_N256.mat');
% load('Pe100eachData_PeFixed_Omega_4100_Nx_1_H11_N256.mat');
eps = [];
h = 10^(-8);
order = 1;
h_con = [];
N = para.HBM.N;
xt = para.xt(:, 6:end);

kn = para.params.func.fc.kn;
xn0 = para.params.func.fc.xn0;
mu = para.params.func.fc.mu;
kt = para.params.func.fc.kt;
w_in = para.params.func.fc.w;

% figure;
% for i = 1:3
%     subplot(3, 3, i);
%     plot(para.t, xt(:, i)), grid on;
%     xlabel('Time');
%     yname = 'x_' + string(i);
%     ylabel(yname);
% end
% 
% for i = 1:3
%     subplot(3, 3, i + 3);
%     plot(para.t, para.Ft(:, i)), grid on;
%     xlabel('Time');
%     yname = 'F_' + string(i);
%     ylabel(yname);
% end
% 
% for i = 1:3
%     subplot(3, 3, i + 6);
%     plot(xt(:, i), para.Ft(:, i)), grid on;
%     xname = 'x_' + string(i);
%     xlabel(xname);
%     yname = 'F_' + string(i);
%     ylabel(yname);
% end

for iH = 1:20
    for iN = 1:12
        
        H = iH;
        ih = iN + 0;
        h = 10^(-ih);
        
        % Numerical Jacobian
        [E, EH] = HBM.fft_matrices(N, H);
        para.params.func.HBM.E = E;
        para.params.func.HBM.EH = EH;
        X = EH * xt;
        xtp = E * X;
        X = X(:);
        % [F, xtF] = fftgx(X, E, EH); % for testing
        JNL_N = finite_diff_jac(@(x) fftgx_f(x, para.params.func), X, h, order);
        
        % Analytiacl Jacobian
        
        [JNL_A, Mft, Gp, dgt] = HBMJACOB_analytical_gf(xt, kn, xn0, mu, kt, w_in, H, para.xp);
        
        % figure;
        % 
        % MftT = [];
        % MftT(:,1) = Mft(3,1,:);
        % subplot(2,2,1)
        % plot(para.t', MftT,'r.'), grid on;
        % xlabel('Time');
        % ylabel('stick condition of x_2');
        % 
        % subplot(2,2,3)
        % MftT(:,1) = Mft(4,1,:);
        % plot(para.t', MftT,'b.'), grid on;
        % xlabel('Time');
        % ylabel('slip condition of x_2');
        % 
        % subplot(2,2,2)
        % dgtT(:,1) = dgt(2,2,:);
        % plot(para.t, dgtT, 'r.'), grid on;
        % xlabel('Time');
        % ylabel('dgt22');
        % 
        % subplot(2,2,4)
        % dgtT(:,1) = dgt(2,3,:);
        % plot(para.t, dgtT, 'b.'), grid on;
        % xlabel('Time');
        % ylabel('dgt23');

        eps(iH, iN) = norm(JNL_A - JNL_N) / norm(JNL_A);
        h_con(iH, iN) = h;

        eps11(iH, iN) = norm(JNL_A(1:(2*H+1), 1:(2*H+1)) - JNL_N(1:(2*H+1), 1:(2*H+1))) / norm(JNL_A(1:(2*H+1), 1:(2*H+1)));
        eps22(iH, iN) = norm(JNL_A((1 + (2*H+1)): 2 * (2*H+1), (1 + (2*H+1)): 2 * (2*H+1)) - JNL_N((1 + (2*H+1)): 2 * (2*H+1), (1 + (2*H+1)): 2 * (2*H+1))) / norm(JNL_A((1 + (2*H+1)): 2 * (2*H+1), (1 + (2*H+1)): 2 * (2*H+1)));
        eps23(iH, iN) = norm(JNL_A((1 + (2*H+1)): 2 * (2*H+1), (1 + 2 * (2*H+1)): 3 * (2*H+1)) - JNL_N((1 + (2*H+1)): 2 * (2*H+1), (1 + 2 * (2*H+1)): 3 * (2*H+1))) / norm(JNL_A((1 + (2*H+1)): 2 * (2*H+1), (1 + 2 * (2*H+1)): 3 * (2*H+1)));
        eps33(iH, iN) = norm(JNL_A((1 + 2 * (2*H+1)): 3 * (2*H+1), (1 + 2 * (2*H+1)): 3 * (2*H+1)) - JNL_N((1 + 2 * (2*H+1)): 3 * (2*H+1), (1 + 2 * (2*H+1)): 3 * (2*H+1))) / norm(JNL_A((1 + 2 * (2*H+1)): 3 * (2*H+1), (1 + 2 * (2*H+1)): 3 * (2*H+1)));
        
    end
end

%% Coulomb friction of dummy fucntion 2 dofs
clear
clc
close all
eps = zeros(11, 10);
h = 10^(-8);
order = 1;
h_con = zeros(11, 10);
N = 256;
H = 1;
A = 0.5;
dt = 2 * pi / N;
t = 0:dt:(2 * pi - dt);
xt = A * fxt(t);
xn = ones(N, 1);
x = [xt(:, 3), xn];

figure; % displacement
plot(t, x), grid on;
legend("x1", "xn");
title('displacement');

kt = 1;
kn = 1;
mu = 0.5;
w = -1.0;
xn0 = 0; % normal pre-displacement

nloop = 2;
[JNL_A, Mft_A, Gp, dgt_A, Ft_A, wt_A] = HBMJACOB_analytical_gf_2dofs(x, kn, xn0, mu, kt, w, H, nloop);
[E, EH] = fft_matrices(N, H);
X = EH * x;
X = X(:);
[JNL_N, Mft_N, dgt_N, Ft_N, wt_N, Ffft] = HBMJACOB_numerical_gf_2dofs(X, kn, xn0, mu, kt, w, H, N, nloop, h, order);

eps = norm(JNL_A - JNL_N) / norm(JNL_A);

figure; % friction forces
plot(t, Ft_A(:, 1)), hold on;
plot(t, mu * Ft_A(:,2), 'k-'), hold on;
plot(t, - mu * Ft_A(:,2), 'k-'), hold on;
legend("T", "mu*Fn", "-mu*Fn");
ylim(1.2 * [min(-mu * Ft_A(:,2)), max(mu * Ft_A(:,2))]);
title('friction forces');
grid on;

figure; % hysteresis cycle
plot(x(:, 1), Ft_A(:, 1));
title('hysteresis cycle');
grid on;

figure; % dgt comparison
subplot(2,2,1)
dgtT(:, 1) = dgt_A(1,1,:);
plot(t, dgtT, 'b-'), hold on;
dgtT(:, 1) = dgt_N(1,1,:);
plot(t, dgtT, 'r--'), hold on;
legend('dgt11_A', 'dgt11_N');
grid on;

subplot(2,2,2)
dgtT(:, 1) = dgt_A(1,2,:);
plot(t, dgtT, 'b-'), hold on;
dgtT(:, 1) = dgt_N(1,2,:);
plot(t, dgtT, 'r--'), hold on;
legend('dgt12_A', 'dgt12_N');
grid on;

subplot(2,2,3)
dgtT(:, 1) = dgt_A(2,1,:);
plot(t, dgtT, 'b-'), hold on;
dgtT(:, 1) = dgt_N(2,1,:);
plot(t, dgtT, 'r--'), hold on;
legend('dgt21_A', 'dgt21_N');
grid on;

subplot(2,2,4)
dgtT(:, 1) = dgt_A(2,2,:);
plot(t, dgtT, 'b-'), hold on;
dgtT(:, 1) = dgt_N(2,2,:);
plot(t, dgtT, 'r--'), hold on;
legend('dgt22_A', 'dgt22_N');
grid on;