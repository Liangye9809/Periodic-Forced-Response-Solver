%% Coulomb friction of dummy fucntion 2 dofs
clear
clc
close all
eps = [];
h = 10^(-2);
order = 1;
h_con = [];
N = 32;
H = 10;
dt = 2 * pi / N;
t = (0:(N-1)) * 2 * pi / N;
t = t';
xn = ones(N, 1);
% xn = - 4 * sin(sin(t)) + 1; % separation to stick
% xt = 2 * exp(cos(t + 1)) - 3; % separation to stick
% xn = 2 * exp(cos(t)) - 0.5; % slip to stick
% xn = 2 * exp(cos(t)) - 0.75; % separation to slip
% xt = 2 * sin(sin(t)); % slip to stick
xt = sin(sin(t)); % pure stick
x = [xt, xn];



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

p.HBM.E = E;
p.HBM.EH = EH;
p.HBM.Nx = 1;
p.HBM.N = N;
p.HBM.H = H;
p.fc.kn = kn;
p.fc.kt = [kt; kt];
p.fc.xn0 = xn0;
p.fc.mu = [mu; mu];
p.fc.w = [w; w];
p.fc.nloop = nloop;

[JNL_A_2, Mft_A_2, JNLt_A, Ft_A_2, wt_A_2] = HBMJACOB_analytical_gf_2dofs_2(x, dx, kn, xn0, mu, kt, w, H, N, nloop);
[JNL_F, JNLt_F, Ft_F, wt_F, Mft_F] = JNL_Analytical(X, p);
function dX = dXinFourier(X, H)
    dX = zeros(size(X));
    for i = 1:H
        dX(2 * i, :) =  i .* X(2 * i + 1, :);
        dX(2 * i + 1, :) =  -i .* X(2 * i, :);
    end

end


%% print
T = t; xplot = x;
for i = 1:nloop - 1
    Tend = T(end);
    T = [T; t + Tend + dt];

    xplot = [xplot; x];
end


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
epsN = norm(JNL_A_2(end-N+1:end,end-N+1:end) - JNL_F(end-N+1:end,end-N+1:end)) / norm(JNL_A_2(end-N+1:end,end-N+1:end));

dTdxt_time_AN = JNLt_A(1:N, 1:2*H+1);
dTdxn_time_AN = JNLt_A(1:N, 2*H+2:end);

% dTdxt_time_F = JNLt_F(1:N, 1:2*H+1);
% dTdxn_time_F = JNLt_F(1:N, 2*(2*H+1)+1:end);

dTdxt_time_F = JNLt_F(N + 1:2 * N, 2*H+1 + 1:2 * (2*H+1));
dTdxn_time_F = JNLt_F(N + 1:2 * N, 2*(2*H+1)+1:end);
% num = 18;
for i = 1:10
    fig = figure; %('PaperOrientation','landscape','PaperUnits','centimeters','PaperPosition', 100 * [0 0 29.7 21], 'PaperSize',[29.7 21] * 100);
    
    fig.WindowState = 'maximized';
    
   
    pause(0.5)
    subplot(2,2,1)
    plot(t, dTdxt_time_F(:, i), 'LineWidth', 2), hold on
    plot(t, dTdxt_time_AN(:, i), 'LineWidth', 2), grid on;
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
    plot(t, dTdxt_time_F(:, i) - dTdxt_time_AN(:, i), 'LineWidth', 2), hold on
    legend('difference');
    grid on

    subplot(2,2,2)
    plot(t, dTdxn_time_F(:, i), 'LineWidth', 2), hold on
    plot(t, dTdxn_time_AN(:, i), 'LineWidth', 2), grid on;
    legend('dTdxn Nu', 'dTdxn An');
    
    title(titlename);

    subplot(2,2,4)
    plot(t, dTdxn_time_F(:, i) - dTdxn_time_AN(:, i), 'LineWidth', 2), hold on
    legend('difference');
    grid on

    % pintname = string(titlename) + '.png';
    % set(gcf,'PaperOrientation','landscape');
    % exportgraphics(gcf, pintname,'ContentType','vector');

end



