%% Coulomb friction of dummy fucntion 2 dofs
clear
clc
% close all
eps = [];
h = 10^(-7);
order = 1;
h_con = [];
N = 2^6;
H = 4;
dt = 2 * pi / N;
t = (0:(N-1)) * 2 * pi / N;
t = t';
xn = ones(N, 1);
% xn = - 4 * sin(sin(t)) + 1; % separation to stick
% xt = 2 * exp(cos(t + 1)) - 3; % separation to stick
% xn = 2 * exp(cos(t)) - 0.5; % slip to stick
% xn = 2 * exp(cos(t)) - 0.75; % separation to slip
% xt = 2 * sin(sin(t)); % slip to stick

% xt = 1.05  * sin(2 .* exp(cos(t))); % tangent case1
% xt = 1.00 * sin(sin(t)) ./ sin(1); % tangent case2
xt = 0.5 * sin(sin(t)) ./ sin(1) + 0.5; % tangent case3 only one side
x = [xt, xn];

% figure; % displacement
% plot(t, x), grid on;
% legend("x1", "xn");
% title('displacement');

kt = 1;
kn = 2;
mu = 0.5;
w =  0;
xn0 = 0; % normal pre-displacement

nloop = 2;

[E, EH] = fft_matrices(N, H);
X = EH * x;
xpr = E * X;
dX = dXinFourier(X, H);
dx = E * dX;
X = X(:);
[Ft_A_2, wt_A_2, Mft_A_2, dxdnt] = gf_2dofs(xpr, kn, xn0, mu, kt, w, nloop);

tic;
        % [JNL_A, Mft_A, Gp, dgt_A, Ft_A, wt_A] = HBMJACOB_analytical_gf_2dofs(x, kn, xn0, mu, kt, w, H, N, nloop);
        [JNL_A_2, JNLt_A] = HBMJACOB_analytical_gf_2dofs_2(dx, kn, mu, kt, H, N, Mft_A_2, dxdnt);
        
toc;
        % eps_a = norm(JNL_A - JNL_A_2) / norm(JNL_A_2);

tic;
        % [JNL_N, Mft_N, dgt_N, Ft_N, wt_N, Ffft] = HBMJACOB_numerical_gf_2dofs(X, kn, xn0, mu, kt, w, H, N, nloop, h, order);
        [JNL_N_2, JNLt_N] = HBMJACOB_numerical_gf_2dofs_2(X, kn, xn0, mu, kt, w, H, N, nloop, h, order);
        % eps_n = norm(JNL_N - JNL_N_2) / norm(JNL_N_2);
toc;

eps_22 = norm(JNL_A_2 - JNL_N_2) / norm(JNL_A_2);

function dX = dXinFourier(X, H)
    dX = zeros(size(X));
    for i = 1:H
        dX(2 * i, :) =  i .* X(2 * i + 1, :);
        dX(2 * i + 1, :) =  -i .* X(2 * i, :);
    end

end
%% calculate G(X+h)
a(:, 1) = Mft_A_2(1, 1, end - N + 1:end);
c(:, 1) = Mft_A_2(3, 1, end - N + 1:end);
Ftemp = fftFt(X, kn, xn0, mu, kt, wt_A_2(end), H, N, nloop);
A_Ft = [a, Ftemp(1:N)];
A_Fn = [c, Ftemp(N + 1:end)];
DX = zeros(2 * (2*H+1), 1);
h = 10^(-7);
for i = 1:size(X, 1)
    Xplus = X;
    Xplus(i) = Xplus(i) + h;
    F_plus_i = fftFt(Xplus, kn, xn0, mu, kt, wt_A_2(end), H, N, nloop);

    A_Ft = [A_Ft, F_plus_i(1:N)];
    A_Fn = [A_Fn, F_plus_i(N + 1:end)];
end
DFt = (A_Ft(:, 3:end) - A_Ft(:, 2)) ./ h;
DFn = (A_Fn(:, 3:end) - A_Fn(:, 2)) ./ h;

% check the xn after cos0 perturbation
XN = X(2*H+2:end);
XNt = E * XN;
figure;
plot(t, xn, 'k-'), hold on;
plot(t, XNt, 'b-'), hold on;
% legend('xn', 'XNt');
XNp = XN;
XNp(1) = XNp(1) + h;
XNpt = E * XNp;
plot(t, XNpt, 'r--'), grid on;
legend('xn', 'XNt', 'XNpt');

dXNt = XNpt - XNt;

Fnt = NormalForces(xn, kn, 0);
FNt = NormalForces(XNt, kn, 0);
FNpt = NormalForces(XNpt, kn, 0);

%%
T = t; xplot = x;
for i = 1:nloop - 1
    Tend = T(end);
    T = [T; t + Tend + dt];

    xplot = [xplot; x];
end

figure;
Mftplot(:, 1) = Mft_A_2(1,1,end-N+1:end);
Mftplot(:, 2) = Mft_A_2(2,1,end-N+1:end);
Mftplot(:, 3) = Mft_A_2(3,1,end-N+1:end);
subplot(2,1,1)
plot(t, Mftplot(:,1), '.'), grid on;
% xlim(T([end-N+1,end]));
title('stick condition');
ylim([-1.2, 1.2]);

subplot(2,1,2)
plot(t, Mftplot(:,3), '.'), grid on;
% xlim(T([end-N+1,end]));
title('contact condition');
ylim([-1.2, 1.2]);

figure; % friction forces
plot(T, Ft_A_2(:, 1), 'LineWidth', 2), hold on;
plot(T, mu * Ft_A_2(:,2), 'k-', 'LineWidth', 2), hold on;
plot(T, - mu * Ft_A_2(:,2), 'k-', 'LineWidth', 2), hold on;
legend("T", "mu*Fn", "-mu*Fn");
ylim(1.2 * [min(-mu * Ft_A_2(:,2)), max(mu * Ft_A_2(:,2))]);
xlim(T([end-N+1,end]));
title('friction forces');
grid on;

figure; % hysteresis cycle
plot(xplot(:, 1), Ft_A_2(:, 1), 'LineWidth', 2);
title('hysteresis cycle');
grid on;


figure;

wp = xplot(:, 1) + mu * kn / kt * max(xplot(:, 2), 0);
wm = xplot(:, 1) - mu * kn / kt * max(xplot(:, 2), 0);
plot(T, wp, 'k--', 'LineWidth', 2), hold on;
plot(T, wm, 'r--', 'LineWidth', 2), hold on;
plot(T, wt_A_2, 'b-', 'LineWidth', 2), hold on;
grid on;
% legend('w+', 'w-', 'wt');

plot(t, x, 'LineWidth', 2), hold on;
legend('w+', 'w-', 'wt', "x1", "xn");
% title('kn = 2, kt = 1, mu = 0.5, xt = 2*sin(sint), xn = 1 / (1 + (cos(t))^2) - 0.5');

figure;

wp = xt + mu * kn / kt * max(xn, 0);
wm = xt - mu * kn / kt * max(xn, 0);
plot(t, wp, 'k--', 'LineWidth', 2), hold on;
plot(t, wm, 'k--', 'LineWidth', 2), hold on;
plot(t, wt_A_2(end-N+1:end), 'b-', 'LineWidth', 2), hold on;
grid on;
legend('w+', 'w-', 'wt');

figure;
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

nexttile
wp = xplot(:, 1) + mu * kn / kt * max(xplot(:, 2), 0);
wm = xplot(:, 1) - mu * kn / kt * max(xplot(:, 2), 0);
plot(T, wp, 'k--', 'LineWidth', 2), hold on;
plot(T, wm, 'r--', 'LineWidth', 2), hold on;
plot(T, wt_A_2, 'b-', 'LineWidth', 2), hold on;
grid on;
% xlim(T([end-2*N+1,end]));
xlim(T([end-N+1,end]));

nexttile
plot(T, Ft_A_2(:, 1), 'LineWidth', 2), hold on;
plot(T, mu * Ft_A_2(:,2), 'k-', 'LineWidth', 2), hold on;
plot(T, - mu * Ft_A_2(:,2), 'k-', 'LineWidth', 2), hold on;
legend("T", "mu*Fn", "-mu*Fn");
ylim(1.2 * [min(-mu * Ft_A_2(:,2)), max(mu * Ft_A_2(:,2))]);
% xlim(T([end-2*N+1,end]));
xlim(T([end-N+1,end]));
grid on;
%% print
figure;
% subplot(3,1,1)
wp = xplot(:, 1) + mu * kn / kt * max(xplot(:, 2), 0);
wm = xplot(:, 1) - mu * kn / kt * max(xplot(:, 2), 0);
plot(T, wp, 'k--', 'LineWidth', 2), hold on;
plot(T, wm, 'r--', 'LineWidth', 2), hold on;
plot(T, wt_A_2, 'b-', 'LineWidth', 2), hold on;
grid on;
plot(t, x, 'LineWidth', 2), hold on;
legend('w+', 'w-', 'wt', "x1", "xn");
title('displacement');

figure;
% subplot(3,1,2)
plot(xplot(:, 1), Ft_A_2(:, 1), 'LineWidth', 2);
title('hysteresis cycle');
grid on;

figure;
% subplot(3,1,3)
plot(T, Ft_A_2(:, 1), 'LineWidth', 2), hold on;
plot(T, mu * Ft_A_2(:,2), 'k-', 'LineWidth', 2), hold on;
plot(T, - mu * Ft_A_2(:,2), 'k-', 'LineWidth', 2), hold on;
legend("T", "mu*Fn", "-mu*Fn");
ylim(1.2 * [min(-mu * Ft_A_2(:,2)), max(mu * Ft_A_2(:,2))]);
title('friction forces');
grid on;

%%
dTdxt_time_NUM = JNLt_N(1:N, 1:2*H+1);
dTdxn_time_NUM = JNLt_N(1:N, 2*H+2:end);

dTdxt_time_AN = JNLt_A(1:N, 1:2*H+1);
dTdxn_time_AN = JNLt_A(1:N, 2*H+2:end);

% num = 18;
for i = 1:1:1
    fig = figure; %('PaperOrientation','landscape','PaperUnits','centimeters','PaperPosition', 100 * [0 0 29.7 21], 'PaperSize',[29.7 21] * 100);
    
    fig.WindowState = 'maximized';
    
    % set(fig,'PaperUnits','centimeters');
    % set(fig,'PaperSize',[29.7 21]);
    % set(fig,'PaperPosition',[0 0 29.7 21]);
    pause(0.5)
    subplot(2,2,1)
    plot(t, dTdxt_time_NUM(:, i), 'LineWidth', 2), hold on
    grid on
    plot(t, dTdxt_time_AN(:, i), 'LineWidth', 2)
    legend('dTdxt Nu', 'dTdxt An');
    if i == 1
        titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'cos0';
    elseif mod(i, 2) == 0
        titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'cos' + string(floor(i/2));
    else
        titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'sin' + string(floor(i/2));
    end
    title(titlename);

    subplot(2,2,3)
    plot(t, dTdxt_time_NUM(:, i) - dTdxt_time_AN(:, i), 'LineWidth', 2), hold on
    legend('difference');
    grid on

    subplot(2,2,2)
    plot(t, dTdxn_time_NUM(:, i), 'LineWidth', 2), hold on
    grid on
    plot(t, dTdxn_time_AN(:, i), 'LineWidth', 2)
    legend('dTdxn Nu', 'dTdxn An');
    
    title(titlename);

    subplot(2,2,4)
    plot(t, dTdxn_time_NUM(:, i) - dTdxn_time_AN(:, i), 'LineWidth', 2), hold on
    legend('difference');
    grid on

    % pintname = string(titlename) + '.png';
    % set(gcf,'PaperOrientation','landscape');
    % exportgraphics(gcf, pintname,'ContentType','vector');

end



