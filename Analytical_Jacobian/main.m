%% Coulomb friction of dummy fucntion 2 dofs
clear
clc
% close all
eps = [];
h = 10^(-4);
order = 1;
h_con = [];
N = 2048;
H = 10;
dt = 2 * pi / N;
t = (0:(N-1)) * 2 * pi / N;
t = t';
% xn = 2*ones(N, 1);
% xn = - 4 * sin(sin(t)) + 1; % separation to stick
% xt = 2 * exp(cos(t + 1)) - 3; % separation to stick
% xn = 2 * exp(cos(t)) - 0.5; % slip to stick
% xn = 2 * exp(cos(t)) - 0.75; % separation to slip
% xt = 2 * sin(sin(t)); % slip to stick
% xt = sin(sin(t)); % pure stick

% plot w for the paper
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
xn = - 4 * sin(sin(t)) + 1; % separation to stick
xt = 2 * exp(cos(t + 1)) - 3; % separation to stick
% simple x
% xn = - 10 * cos(t) + 3; % separation to stick
% xt = - sin(t); % separation to stick

x = [xt, xn];

figure;
plot(t, xt, 'b-'), hold on
plot(t, xn, 'r-');

kt = 1;
kn = 2;
mu = 0.5;
w =  -1;
xn0 = 0; % normal pre-displacement

nloop = 2;

[E, EH] = fft_matrices(N, H);
X = EH * x;
dX = dXinFourier(X, H);
dx = E * dX;
X = [EH * x(:, 1), X];
X = X(:);

[Ft, wt, Mft, dxdnt] = gf_2dofs(x, kn, xn0, mu, kt, w, nloop);
[JNL_A, JNLt_A] = HBMJACOB_analytical_gf_2dofs_2(dx, kn, mu, kt, H, N, Mft, dxdnt);

dxdn = dxdnt(find(dxdnt(:, 1) > 0), 2);
JNL_W = HBMJACOB_analytical_W_gf_2dofs_2(dx, kn, mu, kt, H, N, Mft, dxdn);

[JNL_A_2, JNLt_A_2] = HBMJACOB_analytical_in_time(kn, mu, kt, H, N, Mft, dxdn);

function dX = dXinFourier(X, H)
    dX = zeros(size(X));
    for i = 1:H
        dX(2 * i, :) =  i .* X(2 * i + 1, :);
        dX(2 * i + 1, :) =  -i .* X(2 * i, :);
    end

end

% eps1 = norm(JNL_W - JNL_A) / norm(JNL_W)
% eps2 = norm(JNL_W - JNL_A_2) / norm(JNL_W)
% eps3 = norm(JNL_A_2 - JNL_A) / norm(JNL_A_2)

% print
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
set(gca,'TickLabelInterpreter','latex')
%%
figure;
% subplot(3,1,2)
plot(xplot(:, 1), Ft(:, 1), 'LineWidth', 2);
title('hysteresis cycle');
grid on;

figure;
% subplot(3,1,3)
plot(T, Ft(:, 1), 'LineWidth', 2), hold on;
plot(T, mu * Ft(:,2), 'k-', 'LineWidth', 2), hold on;
plot(T, - mu * Ft(:,2), 'k-', 'LineWidth', 2), hold on;
legend("T", "mu*Fn", "-mu*Fn");
ylim(1.2 * [min(-mu * Ft(:,2)), max(mu * Ft(:,2))]);
title('friction forces');
grid on;

%%
% epsN = norm(JNLt_A_2(end-N+1:end,end-N+1:end) - JNL_F(end-N+1:end,end-N+1:end)) / norm(JNL_A_2(end-N+1:end,end-N+1:end));

dTdxt_time_AN = JNLt_A(1:N, 1:2*H+1);
dTdxn_time_AN = JNLt_A(1:N, 2*H+2:end);

dTdxt_time_AN_2 = JNLt_A_2(1:N, 1:2*H+1);
dTdxn_time_AN_2 = JNLt_A_2(1:N, 2*H+2:end);

% num = 18;
for i = 1:3
    fig = figure; %('PaperOrientation','landscape','PaperUnits','centimeters','PaperPosition', 100 * [0 0 29.7 21], 'PaperSize',[29.7 21] * 100);
    
    fig.WindowState = 'maximized';

    subplot(2,2,1)
    plot(t, dTdxt_time_AN_2(:, i), 'LineWidth', 2), hold on
    plot(t, dTdxt_time_AN(:, i), 'LineWidth', 2), grid on;
    legend('dTdxt 2', 'dTdxt 1');
    if i == 1
        titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'cos0';
    elseif mod(i, 2) == 0
        titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'cos' + string(floor(i/2));
    else
        titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'sin' + string(floor(i/2));
    end
    title(titlename);

    subplot(2,2,3)
    plot(t, dTdxt_time_AN_2(:, i) - dTdxt_time_AN(:, i), 'LineWidth', 2), hold on
    legend('difference');
    grid on

    subplot(2,2,2)
    plot(t, dTdxt_time_AN_2(:, i), 'LineWidth', 2), hold on
    plot(t, dTdxn_time_AN(:, i), 'LineWidth', 2), grid on;
    legend('dTdxn 2', 'dTdxn 1');
    
    title(titlename);

    subplot(2,2,4)
    plot(t, dTdxt_time_AN_2(:, i) - dTdxn_time_AN(:, i), 'LineWidth', 2), hold on
    legend('difference');
    grid on

    % pintname = string(titlename) + '.png';
    % set(gcf,'PaperOrientation','landscape');
    % exportgraphics(gcf, pintname,'ContentType','vector');

end


%%
x = 4;
y1 = 0.91;
y2 = 1.09;

quiver(x,y1,0,y2-y1,...
    0,...
    'MaxHeadSize',0.3,...
    'LineWidth',1.5)