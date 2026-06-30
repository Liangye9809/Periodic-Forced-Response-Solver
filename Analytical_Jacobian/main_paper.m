%% Coulomb friction of dummy fucntion 2 dofs
clear
clc
% close all
eps = [];
h = 10^(-7);
order = 1;
h_con = [];
N = 256;
H = 1;
dt = 2 * pi / N;
t = (0:(N-1)) * 2 * pi / N;
t = t';
% xn = ones(N, 1);
% xn = - 4 * sin(sin(t)) + 1; % separation to stick
% xt = 2 * exp(cos(t + 1)) - 3; % separation to stick
% xn = 2 * exp(cos(t)) - 0.5; % slip to stick
% xn = 2 * exp(cos(t)) - 0.75; % separation to slip
% xt = 2 * sin(sin(t)); % slip to stick

% xt = 1.05  * sin(2 .* exp(cos(t))); % tangent case1
% xt = 1.00 * sin(sin(t)) ./ sin(1); % tangent case2
% xt = 0.5 * sin(sin(t)) ./ sin(1) + 0.5; % tangent case3 only one side

% plot w for the paper (kt = 1, kn = 2, mu = 0.5)
% pure stick
% xn = 2*ones(N, 1); % pure stick
% xt = sin(sin(t)); % pure stick
% simple x
% xn = 2*ones(N, 1); % pure stick
% xt = sin(t); % pure stick

% slip to stick
% xn = 2 * exp(cos(t)) - 0.5; % slip to stick
% xt = 2 * sin(sin(t)); % slip to stick
% simple x
% xn = 2.5 * cos(t) + 3; % slip to stick
% xt = 3 * sin(t); % slip to stick

% gap to stick
% xn = - 4 * sin(sin(t)) + 1; % separation to stick
% xt = 2 * exp(cos(t + 1)) - 3; % separation to stick
% simple x
% xn = - 10 * cos(t) + 3; % separation to stick
% xt = - sin(t); % separation to stick
xn = - 1 * cos(t) + 0.5; % separation to stick % (kt = 1, kn = 500, mu = 0.8)
xt = 2 * cos(t + pi/6); % separation to stick
x = [xt, xn];

figure; % displacement
plot(t, x, 'LineWidth', 2), grid on;
legend("x1", "xn");
title('displacement');
xticks([0:pi/2:2*pi]);
xticklabels({'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
set(gca,'TickLabelInterpreter','latex');

kt = 1;
kn = 500;
mu = 0.8;
w =  0;
xn0 = 0; % normal pre-displacement

nloop = 2;

[E, EH] = fft_matrices(N, H);
X = EH * x;
xpr = E * X;
dX = dXinFourier(X, H);
dx = E * dX;
X = X(:);
[Ft, wt, Mft, dxdnt] = gf_2dofs(xpr, kn, xn0, mu, kt, w, nloop);

tic;
        [JNL_A, JNLt_A] = HBMJACOB_analytical_gf_2dofs_2(dx, kn, mu, kt, H, N, Mft, dxdnt);
        
toc;
        % eps_a = norm(JNL_A - JNL_A_2) / norm(JNL_A_2);

tic;
        [JNL_N, JNLt_N] = HBMJACOB_numerical_gf_2dofs_2(X, kn, xn0, mu, kt, w, H, N, nloop, h, order);
        % eps_n = norm(JNL_N - JNL_N_2) / norm(JNL_N_2);
toc;


eps_22 = norm(JNL_A - JNL_N) / norm(JNL_A);
dxdn = dxdnt(find(dxdnt(:,1) == 1),2)
function dX = dXinFourier(X, H)
    dX = zeros(size(X));
    for i = 1:H
        dX(2 * i, :) =  i .* X(2 * i + 1, :);
        dX(2 * i + 1, :) =  -i .* X(2 * i, :);
    end

end
%% calculate G(X+h)
% a(:, 1) = Mft(1, 1, end - N + 1:end);
% c(:, 1) = Mft(3, 1, end - N + 1:end);
% Ftemp = fftFt(X, kn, xn0, mu, kt, wt(end), H, N, nloop);
% A_Ft = [a, Ftemp(1:N)];
% A_Fn = [c, Ftemp(N + 1:end)];
% DX = zeros(2 * (2*H+1), 1);
% h = 10^(-7);
% for i = 1:size(X, 1)
%     Xplus = X;
%     Xplus(i) = Xplus(i) + h;
%     F_plus_i = fftFt(Xplus, kn, xn0, mu, kt, wt(end), H, N, nloop);
% 
%     A_Ft = [A_Ft, F_plus_i(1:N)];
%     A_Fn = [A_Fn, F_plus_i(N + 1:end)];
% end
% DFt = (A_Ft(:, 3:end) - A_Ft(:, 2)) ./ h;
% DFn = (A_Fn(:, 3:end) - A_Fn(:, 2)) ./ h;
% 
% % check the xn after cos0 perturbation
% XN = X(2*H+2:end);
% XNt = E * XN;
% figure;
% plot(t, xn, 'k-'), hold on;
% plot(t, XNt, 'b-'), hold on;
% % legend('xn', 'XNt');
% XNp = XN;
% XNp(1) = XNp(1) + h;
% XNpt = E * XNp;
% plot(t, XNpt, 'r--'), grid on;
% legend('xn', 'XNt', 'XNpt');
% 
% dXNt = XNpt - XNt;
% 
% Fnt = NormalForces(xn, kn, 0);
% FNt = NormalForces(XNt, kn, 0);
% FNpt = NormalForces(XNpt, kn, 0);

%%
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultAxesFontSize',16)
set(0, 'DefaultFigurePosition', [500, 500, 600, 450]);
set(0, 'DefaultFigureColor', 'w');

T = t; xplot = x;
for i = 1:nloop - 1
    Tend = T(end);
    T = [T; t + Tend + dt];

    xplot = [xplot; x];
end

figure;
Mftplot(:, 1) = Mft(1,1,end-N+1:end);
Mftplot(:, 2) = Mft(2,1,end-N+1:end);
Mftplot(:, 3) = Mft(3,1,end-N+1:end);
subplot(3,1,1)
plot(t, Mftplot(:,1), '.'), grid on;
% xlim(T([end-N+1,end]));
title('stick condition');
ylim([-1.2, 1.2]);

subplot(3,1,2)
plot(t, Mftplot(:,2), '.'), grid on;
% xlim(T([end-N+1,end]));
title('slip condition');
ylim([-1.2, 1.2]);

subplot(3,1,3)
plot(t, Mftplot(:,3), '.'), grid on;
% xlim(T([end-N+1,end]));
title('contact condition');
ylim([-1.2, 1.2]);

figure; % friction forces
plot(t, Ft(end-N+1:end, 1), 'LineWidth', 2), hold on;
plot(t, mu * Ft(end-N+1:end,2), 'k-', 'LineWidth', 2), hold on;
plot(t, - mu * Ft(end-N+1:end,2), 'k-', 'LineWidth', 2), hold on;
legend("T", "mu*Fn", "-mu*Fn");
% ylim(1.2 * [min(-mu * Ft(:,2)), max(mu * Ft(:,2))]);
% xlim(T([end-N+1,end]));
title('friction forces');
grid on;
xticks([0:pi/2:2*pi]);
xticklabels({'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
set(gca,'TickLabelInterpreter','latex');

figure; % hysteresis cycle
plot(xplot(:, 1), Ft(:, 1), 'LineWidth', 2);
title('hysteresis cycle');
grid on;


figure; % ('Units','centimeters','Position',[2, 2, 6, 4.5]);
% subplot(3,1,1)
wp = xplot(:, 1) + mu * kn / kt * max(xplot(:, 2), 0);
wm = xplot(:, 1) - mu * kn / kt * max(xplot(:, 2), 0);
w_max = min(wp);
w_min = max(wm);
plot(T, wp, 'k--', 'LineWidth', 2, 'DisplayName', '$w^+$'), hold on;
plot(T, wm, 'r--', 'LineWidth', 2, 'DisplayName', '$w^-$'), hold on;
plot(T, wt, 'b-', 'LineWidth', 2, 'DisplayName', '$w$'), hold on;
grid on;
% plot(t, x, 'LineWidth', 2), hold on;
% plot(T, [x(:,1); x(:,1)], 'LineWidth', 2, 'DisplayName', 'x1'), hold on;
% plot(T, [x(:,2); x(:,2)], 'LineWidth', 2, 'DisplayName', 'xn'), hold on;
% plot(T, w_max * ones(size(T)), 'b-', 'LineWidth', 2, 'HandleVisibility','off');
% plot(T, w_min * ones(size(T)), 'b-', 'LineWidth', 2, 'HandleVisibility','off');
legend show;
% annotation('doublearrow', [0.25, 0.25], [(w_min + 5)/10, (w_max + 5)/10]);
% quiver(5, 0, 0, w_max, 0, 'MaxHeadSize', 0.5, 'LineWidth', 1.5, 'Color', 'b', 'HandleVisibility','off');
% quiver(5, 0, 0, w_min, 0, 'MaxHeadSize', 0.5, 'LineWidth', 1.5, 'Color', 'b', 'HandleVisibility','off');
% ylim([-3, 3]);
xlim([0, 4*pi]);
% text(5.5,0,'$w_0$','Interpreter','latex', 'Color','b', 'FontSize', 16);
% title('displacement');
xlabel('Time');
ylabel('Displacement');
xticks([0:pi/2:4*pi]);
xticklabels({'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$',...
    '$\frac{5\pi}{2}$','$3\pi$','$\frac{7\pi}{2}$','$4\pi$'})
set(gca,'TickLabelInterpreter','latex');



figure;
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

nexttile
wp = xplot(:, 1) + mu * kn / kt * max(xplot(:, 2), 0);
wm = xplot(:, 1) - mu * kn / kt * max(xplot(:, 2), 0);
plot(T, wp, 'k--', 'LineWidth', 2), hold on;
plot(T, wm, 'r--', 'LineWidth', 2), hold on;
plot(T, wt, 'b-', 'LineWidth', 2), hold on;
grid on;
% xlim(T([end-2*N+1,end]));
xlim(T([end-N+1,end]));

nexttile
plot(T, Ft(:, 1), 'LineWidth', 2), hold on;
plot(T, mu * Ft(:,2), 'k-', 'LineWidth', 2), hold on;
plot(T, - mu * Ft(:,2), 'k-', 'LineWidth', 2), hold on;
legend("T", "mu*Fn", "-mu*Fn");
ylim(1.2 * [min(-mu * Ft(:,2)), max(mu * Ft(:,2))]);
% xlim(T([end-2*N+1,end]));
xlim(T([end-N+1,end]));
grid on;



%%
dTdxt_time_NUM = JNLt_N(1:N, 1:2*H+1);
dTdxn_time_NUM = JNLt_N(1:N, 2*H+2:end);

dTdxt_time_AN = JNLt_A(1:N, 1:2*H+1);
dTdxn_time_AN = JNLt_A(1:N, 2*H+2:end);

% num = 18;
for i = 1:3
    fig = figure; %('PaperOrientation','landscape','PaperUnits','centimeters','PaperPosition', 100 * [0 0 29.7 21], 'PaperSize',[29.7 21] * 100);
    
    fig.WindowState = 'maximized';
    
    % set(fig,'PaperUnits','centimeters');
    % set(fig,'PaperSize',[29.7 21]);
    % set(fig,'PaperPosition',[0 0 29.7 21]);
    
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

%% compare 3 Jacobians in time
clear
J_AN = load('P_gap_to_stick_Analytical_and_Numerical.mat');
J_F  = load('P_gap_to_stick_Fixed_Numerical.mat');

JNLt_N = J_AN.P.JNLt_N;
JNLt_A = J_AN.P.JNLt_A;
JNLt_F = J_F.P.JNLt_N;
N = J_F.P.N;
H = J_F.P.H;
t = (0:(N-1)) * 2 * pi / N;
t = t';
dTdxt_time_NUM = JNLt_N(1:N, 1:2*H+1);   dNdxt_time_NUM = JNLt_N(N + 1:end, 1:2*H+1);
dTdxn_time_NUM = JNLt_N(1:N, 2*H+2:end); dNdxn_time_NUM = JNLt_N(N + 1:end, 2*H+2:end);

dTdxt_time_AN = JNLt_A(1:N, 1:2*H+1);   dNdxt_time_AN = JNLt_A(N + 1:end, 1:2*H+1);
dTdxn_time_AN = JNLt_A(1:N, 2*H+2:end); dNdxn_time_AN = JNLt_A(N + 1:end, 2*H+2:end);

dTdxt_time_F = JNLt_F(1:N, 1:2*H+1);   dNdxt_time_F = JNLt_F(N + 1:end, 1:2*H+1);
dTdxn_time_F = JNLt_F(1:N, 2*H+2:end); dNdxn_time_F = JNLt_F(N + 1:end, 2*H+2:end);

for i = 1:3
    fig = figure; %('PaperOrientation','landscape','PaperUnits','centimeters','PaperPosition', 100 * [0 0 29.7 21], 'PaperSize',[29.7 21] * 100);
    
    fig.WindowState = 'maximized';
    
    % set(fig,'PaperUnits','centimeters');
    % set(fig,'PaperSize',[29.7 21]);
    % set(fig,'PaperPosition',[0 0 29.7 21]);
    
    subplot(2,2,1)
    plot(t, dTdxt_time_AN(:, i),'k-', 'LineWidth', 2, 'DisplayName', 'Analytical'), hold on   
    plot(t, dTdxt_time_NUM(:, i),'g-.', 'LineWidth', 2, 'DisplayName', 'Numerical'), grid on
    plot(t, dTdxt_time_F(:, i),'r--', 'LineWidth', 2, 'DisplayName', 'Corrected Numerical')
    legend show;
    if i == 1
        titlename = 'c}^0}$';
        ylim([-1, 1]);
    elseif mod(i, 2) == 0
        titlename = 'c}^' + string(floor(i/2)) + '}$';
    else
        titlename = 's}^' + string(floor(i/2)) + '}$';
    end
    % titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' + string(titlename);
    titlename = '$\frac{\partial f_t(t)}{\partial \mathbf{X}_{t' + string(titlename);
    title(titlename);
    xticks([0:pi/2:2*pi]);
    xticklabels({'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,2)
    plot(t, dTdxn_time_AN(:, i),'k-',  'LineWidth', 2, 'DisplayName', 'Analytical'), hold on    
    plot(t, dTdxn_time_NUM(:, i),'g-.', 'LineWidth', 2, 'DisplayName', 'Numerical'), grid on
    plot(t, dTdxn_time_F(:, i),'r--',  'LineWidth', 2, 'DisplayName', 'Corrected Numerical')
    legend show;
    if i == 1
        titlename = 'c}^0}$';
    elseif mod(i, 2) == 0
        titlename = 'c}^' + string(floor(i/2)) + '}$';
    else
        titlename = 's}^' + string(floor(i/2)) + '}$';
    end
    % titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' + string(titlename);
    titlename = '$\frac{\partial f_t(t)}{\partial \mathbf{X}_{n' + string(titlename);
    title(titlename);
    xticks([0:pi/2:2*pi]);
    xticklabels({'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,3)
    plot(t, dNdxt_time_AN(:, i),'k-',  'LineWidth', 2, 'DisplayName', 'Analytical'), hold on    
    plot(t, dNdxt_time_NUM(:, i),'g-.', 'LineWidth', 2, 'DisplayName', 'Numerical'), grid on
    plot(t, dNdxt_time_F(:, i),'r--',  'LineWidth', 2, 'DisplayName', 'Corrected Numerical')
    legend show;
    if i == 1
        titlename = 'c}^0}$';
    elseif mod(i, 2) == 0
        titlename = 'c}^' + string(floor(i/2)) + '}$';
    else
        titlename = 's}^' + string(floor(i/2)) + '}$';
    end
    % titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' + string(titlename);
    titlename = '$\frac{\partial f_n(t)}{\partial \mathbf{X}_{t' + string(titlename);
    title(titlename);
    xticks([0:pi/2:2*pi]);
    xticklabels({'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,4)
    plot(t, dNdxn_time_AN(:, i),'k-',  'LineWidth', 2, 'DisplayName', 'Analytical'), hold on    
    plot(t, dNdxn_time_NUM(:, i),'g-.', 'LineWidth', 2, 'DisplayName', 'Numerical'), grid on
    plot(t, dNdxn_time_F(:, i),'r--',  'LineWidth', 2, 'DisplayName', 'Corrected Numerical')
    legend show;
    if i == 1
        titlename = 'c}^0}$';
    elseif mod(i, 2) == 0
        titlename = 'c}^' + string(floor(i/2)) + '}$';
    else
        titlename = 's}^' + string(floor(i/2)) + '}$';
    end
    % titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' + string(titlename);
    titlename = '$\frac{\partial f_n(t)}{\partial \mathbf{X}_{n' + string(titlename);
    title(titlename);
    xticks([0:pi/2:2*pi]);
    xticklabels({'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    set(gca,'TickLabelInterpreter','latex');

    % pintname = string(titlename) + '.png';
    % set(gcf,'PaperOrientation','landscape');
    % exportgraphics(gcf, pintname,'ContentType','vector');

end
