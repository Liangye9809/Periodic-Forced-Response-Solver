%% plot Amplitude vs omega


[E, EH] = HBM.fft_matrices(N, H);
for i = 1:3*Nx + Na
    x_contNx(:,:,i) = x_cont((2*H+1)*(i-1)+1:(2*H+1)*i,:);
end
Adof = omega_cont';
for Ndof = 1:3*Nx + Na
    A = E * x_contNx(:,:,Ndof);
    Amax = max(abs(A));
    Adof(:, Ndof + 1) = Amax';
end
Adof = [Adof, k_cont'];

% ind_gap_stick = find(gap_cont == 1 & (slipP_cont + slipM_cont) == 0);
% ind_gap = find(gap_cont == 1);
% ind_slip = find(slipP_cont == 1);

figure;
% yyaxis left
plot(Adof(:, 1), Adof(:, 5), 'r.', 'LineWidth', 2), hold on;
plot(Adof(:, 1), Adof(:, 3), 'b.', 'LineWidth', 2), grid on;
xlabel('Omega');
legend('xn', 'xt');
% ylim([0.1, 100]);
titlename = 'mu = ' + string(mu(1)) + ', kt = ' + string(kt(1)) + ', kn = ' + string(kn) + ', max(xt) = ' + string(max(Adof(:, 3)));
title(titlename);
% yyaxis right
% plot(Adof(:, 1), gap_cont', 'LineWidth', 2, 'LineStyle', '-', 'Color', 'g'), hold on;
% plot(Adof(:, 1), slipP_cont', 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'), hold on;
% stem(Adof(:, 1), k_cont'), hold on;
% ylim([0, 1.1]);
% xlabel('Omega');
% legend('xt,damper 0.2', 'xn,damper 0.4', 'gap appears area', 'slip appears area');

% filename = 'data/Analytical Petrov System 1/ky = ' + string(kn) + ', g = ' + string(xn0);
% save(filename, 'Adof');

% filename = 'data\Analytical Petrov System 1\ky = ' + string(kn) + ', g = ' + string(xn0);
% save(filename, 'Adof');

%%
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultAxesFontSize',16)
set(0, 'DefaultFigurePosition', [500, 500, 600, 450]);
set(0, 'DefaultFigureColor', 'w');




%%
% A_p = load('data/Analytical Petrov System 3/para_H10_ds0.05_N512_stopped_point.mat');

i_plot = 115; % point before second 11.6743
x_poss = x_cont(:, i_plot);
% x_poss = x_cont(:, end);
% omega_poss = omega_cont(i_plot);
% omega_poss = omega_cont(end);
% omega_poss = 11.48;
omega_poss = 11.96;
params.func.fc.w = w_cont(:, i_plot);
% params.func.fc.w = w_cont(:, end);
for i = 1:Na + 3 * Nx
    r1 = (2 * H + 1) * (i - 1) + 1;
    r2 = (2 * H + 1) * i;
    X(:,i) = x_poss(r1:r2); % reorder in dofs in column
end
xt_poss = params.func.HBM.E * X;
xct_poss = xt_poss(:, Na + 1:end);

[FUN_poss, w_poss, ~, flag_poss] = HBMFUNC(x_poss, xct_poss + xp', omega_poss, params.func);

% if FUN(x) ~= 0, calculate corresponse value
if norm(FUN_poss) > params.Newton.epsf
    warning('norm(FUN_poss) > params.Newton.epsf');
    params.cont.ds = 0;
    params.cont.omega_0 = omega_poss;
    params.cont.step = 100001;
    params.cont.x0 = x_poss;
    [x_poss, omega_poss, ~, ~, w_poss, ~, FlagState_poss] = cont_step(@HBMFUNC, @HBMJACOB, @HBMJOmega, params);
    params.func.fc.w = w_poss;

    for i = 1:Na + 3 * Nx
        r1 = (2 * H + 1) * (i - 1) + 1;
        r2 = (2 * H + 1) * i;
        X(:,i) = x_poss(r1:r2); % reorder in dofs in column
    end
    xt_poss = params.func.HBM.E * X; 
    xct_poss = xt_poss(:, Na + 1:end);
end


[Ft_poss, w_poss, flag_poss] = g(xct_poss + xp', kn, xn0, mu, kt, params.func.fc.w, nloop); 
JNL_poss = JNL_Analytical(xct_poss, flag_poss(:, :, end - N + 1:end), H, N, kt, kn, mu);
T = 2 * pi / omega_poss;
dt = T / N;
t_poss = [0:dt:(T - dt)]';

para.t = t_poss;
para.xt = xt_poss;
para.Ft = Ft_poss;
para.omega = omega_poss;
para.xp = xp;
para.gxp = gxp;
para.params = params;
para.X = X;
para.flag = flag_poss;
% solution
para.x_poss = x_poss;
para.x_cont = x_cont;
para.k_cont = k_cont;
para.omega_cont = omega_cont;
para.slipM_cont = slipM_cont;
para.slipP_cont = slipP_cont;
para.stick_cont = stick_cont;
para.gap_cont = gap_cont;
para.JNL_poss = JNL_poss;


figure; % displacement
subplot(2,2,1)
yyaxis left
plot(t_poss, xt_poss(:, 2) + xp(1), 'b-', 'LineWidth', 2, 'DisplayName', 'xt'), grid on;
% figure; % friction
yyaxis right
plot(t_poss, Ft_poss(end - N + 1:end, 1), 'r-', 'LineWidth', 2, 'DisplayName', 'T'), grid on;
legend('show');
titlename = 'Omega = ' + string(omega_poss) + ', Amplitude = ' + string(max(abs(xct_poss(:, 1))));
title(titlename);

subplot(2,2,2)
yyaxis left
plot(t_poss, xt_poss(:, 4) + xp(3), 'b-', 'LineWidth', 2, 'DisplayName', 'xn'), grid on;
% figure; % friction
yyaxis right
plot(t_poss, Ft_poss(end - N + 1:end, 3), 'r-', 'LineWidth', 2, 'DisplayName', 'Fn'), grid on;
legend('show');
titlename = 'Omega = ' + string(omega_poss) + ', Amplitude = ' + string(max(abs(xct_poss(:, 3))));
title(titlename);

subplot(2,2,3)
plot(xt_poss(:, 2) + xp(1), Ft_poss(end - N + 1:end, 1), 'b-', 'LineWidth', 2), grid on;
% legend('show');
xlabel('xt');
ylabel('T');

subplot(2,2,4)
plot(t_poss, mu(1) * Ft_poss(end - N + 1:end, 3), 'k-', 'LineWidth', 2, 'DisplayName', 'mu*Fn'), hold on;
plot(t_poss, -mu(1) * Ft_poss(end - N + 1:end, 3), 'k-', 'LineWidth', 2, 'DisplayName', '-mu*Fn'), grid on;
plot(t_poss, Ft_poss(end - N + 1:end, 1), 'b-', 'LineWidth', 2, 'DisplayName', 'T');
legend('show');
xlabel('t');
ylabel('F');
% savename = 'data/Analytical Petrov System 1/ky = 120, g = 10/Omega = ' + string(omega_poss) + ', Amplitude = ' + string(Adof(i_plot, 5)) + '.mat';
% save(savename, 'para');

figure;
a(:, 1) = flag_poss(1, 1, end - N + 1:end);
b(:, 1) = flag_poss(2, 1, end - N + 1:end);
plot(t_poss, a, 'b.'), grid on;
ylim([-1.2, 2.2]);

% figure;
% missing = kt(1) * (xct_poss(475, 1) - xct_poss(474, 1)) / (xct_poss(475, 3) - xct_poss(474, 3))
% mu_min = kt(1) * (xct_poss(347, 1) - w_poss(1, 1, N + 346)) / Ft_poss(347 + N, 3)

%%
[E, EH] = HBM.fft_matrices(N, H);
JNLt_A_11 = E * JNL_poss(1:2 * H + 1, 1:2 * H + 1);
JNLt_A_13 = E * JNL_poss(1:2 * H + 1, 2 * (2 * H + 1) + 1:end);
JNLt_A_22 = E * JNL_poss((2 * H + 1) + 1:2 * (2 * H + 1), (2 * H + 1) + 1:2 * (2 * H + 1));
JNLt_A_23 = E * JNL_poss((2 * H + 1) + 1:2 * (2 * H + 1), 2 * (2 * H + 1) + 1:end);

dt = T / N;
t = [0:dt:(T - dt)]';

for i = 1:3
    fig = figure; %('PaperOrientation','landscape','PaperUnits','centimeters','PaperPosition', 100 * [0 0 29.7 21], 'PaperSize',[29.7 21] * 100);
    
    fig.WindowState = 'maximized';
    if i == 1
        titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'cos0';
    elseif mod(i, 2) == 0
        titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'cos' + string(floor(i/2));
    else
        titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'sin' + string(floor(i/2));
    end
    title(titlename);

    subplot(2,2,1)
    plot(t, JNLt_A_11(:, i), 'LineWidth', 2), hold on
    legend('dTdxt 11');
    title(titlename);
    

    subplot(2,2,2)
    plot(t, JNLt_A_13(:, i), 'LineWidth', 2), hold on
    legend('dTdxn 13');
    title(titlename);
    

    subplot(2,2,3)
    plot(t, JNLt_A_22(:, i), 'LineWidth', 2), hold on
    legend('dTdxt 22');
    title(titlename);
    

    subplot(2,2,4)
    plot(t, JNLt_A_23(:, i), 'LineWidth', 2), hold on
    legend('dTdxn 23');
    title(titlename);
    
    % pintname = string(titlename) + '.png';
    % set(gcf,'PaperOrientation','landscape');
    % exportgraphics(gcf, pintname,'ContentType','vector');

end

%% correct Numerical
clear
% A_p = load('data/Analytical Petrov System 3/kt = 30, mu = 8, k = 40, f = 100sin(t)/para_H5_ds0.05_N256_stopped_point.mat');
% A_p = load('data/Analytical Petrov System 3/kt = 30, mu = 8, k = 40, f = 100sin(t)/para_H5_ds0.05_N128_stopped_point.mat');
A_p = load('data/Analytical Petrov System 3/kt = 30, mu = 8, k = 40, f = 100sin(t)/para_H5_ds0.09_N128_stopped_point_backward_mu10.mat');

% i_plot = 707; % point before second 11.6743
x_poss = A_p.para.X(:);
% omega_poss = omega_cont(i_plot);
% omega_poss = 11.6743;
omega_poss = A_p.para.omega;
params = A_p.para.params;
xp = A_p.para.xp;
Na = A_p.para.params.func.HBM.Na;
Nx = A_p.para.params.func.HBM.Nx;
N = A_p.para.params.func.HBM.N;
H = A_p.para.params.func.HBM.H;
kn = A_p.para.params.func.fc.kn;
xn0 = A_p.para.params.func.fc.xn0;
mu = A_p.para.params.func.fc.mu;
kt = A_p.para.params.func.fc.kt;
nloop = A_p.para.params.func.fc.nloop;
for i = 1:Na + 3 * Nx
    r1 = (2 * H + 1) * (i - 1) + 1;
    r2 = (2 * H + 1) * i;
    X(:,i) = x_poss(r1:r2); % reorder in dofs in column
end
xt_poss = params.func.HBM.E * X;
xct_poss = xt_poss(:, Na + 1:end);


[FUN_poss, w_poss, ~, flag_poss] = HBMFUNC(x_poss, xct_poss, omega_poss, params.func);

% if FUN(x) ~= 0, calculate corresponse value
if norm(FUN_poss) > params.Newton.epsf
    warning('norm(FUN_poss) > params.Newton.epsf');
    params.cont.ds = 0;
    params.cont.omega_0 = omega_poss;
    params.cont.step = 100001;
    params.cont.x0 = x_poss;
    [x_poss, omega_poss, ~, ~, w_poss, ~, ~] = cont_step(@HBMFUNC, @HBMJACOB, @HBMJOmega, params);
    params.func.fc.w = w_poss;

    for i = 1:Na + 3 * Nx
        r1 = (2 * H + 1) * (i - 1) + 1;
        r2 = (2 * H + 1) * i;
        X(:,i) = x_poss(r1:r2); % reorder in dofs in column
    end
    xt_poss = params.func.HBM.E * X; 
    xct_poss = xt_poss(:, Na + 1:end);
end


[Ft_poss, w_poss, flag_poss] = g(xct_poss + xp', kn, xn0, mu, kt, params.func.fc.w, nloop); 


JNL_poss = finite_diff_jac(@(x) fftgx(x, xct_poss, params.func), x_poss((2 * H + 1) * Na + 1:end));

T = 2 * pi / omega_poss;
dt = T / N;
t_poss = [0:dt:(T - dt)]';

para.t = t_poss;
para.xt = xt_poss;
para.Ft = Ft_poss;
para.omega = omega_poss;
para.xp = xp;

para.params = params;
para.X = X;
para.flag = flag_poss;
% solution


para.JNL_poss = JNL_poss;


figure; % displacement
subplot(2,2,1)
yyaxis left
plot(t_poss, xt_poss(:, 2) + xp(1), 'b-', 'LineWidth', 2, 'DisplayName', 'xt'), grid on;
% figure; % friction
yyaxis right
plot(t_poss, Ft_poss(end - N + 1:end, 1), 'r-', 'LineWidth', 2, 'DisplayName', 'T'), grid on;
legend('show');
titlename = 'Omega = ' + string(omega_poss) + ', Amplitude = ' + string(max(abs(xct_poss(:, 1))));
title(titlename);

subplot(2,2,2)
yyaxis left
plot(t_poss, xt_poss(:, 4) + xp(3), 'b-', 'LineWidth', 2, 'DisplayName', 'xn'), grid on;
% figure; % friction
yyaxis right
plot(t_poss, Ft_poss(end - N + 1:end, 3), 'r-', 'LineWidth', 2, 'DisplayName', 'Fn'), grid on;
legend('show');
titlename = 'Omega = ' + string(omega_poss) + ', Amplitude = ' + string(max(abs(xct_poss(:, 3))));
title(titlename);

subplot(2,2,3)
plot(xt_poss(:, 2) + xp(1), Ft_poss(end - N + 1:end, 1), 'b-', 'LineWidth', 2), grid on;
% legend('show');
xlabel('xt');
ylabel('T');

subplot(2,2,4)
plot(t_poss, mu(1) * Ft_poss(end - N + 1:end, 3), 'k-', 'LineWidth', 2, 'DisplayName', 'mu*Fn'), hold on;
plot(t_poss, -mu(1) * Ft_poss(end - N + 1:end, 3), 'k-', 'LineWidth', 2, 'DisplayName', '-mu*Fn'), grid on;
plot(t_poss, Ft_poss(end - N + 1:end, 1), 'b-', 'LineWidth', 2, 'DisplayName', 'T');
legend('show');
xlabel('t');
ylabel('F');
% savename = 'data/Analytical Petrov System 1/ky = 120, g = 10/Omega = ' + string(omega_poss) + ', Amplitude = ' + string(Adof(i_plot, 5)) + '.mat';
% save(savename, 'para');

figure;
a(:, 1) = flag_poss(1, 1, end - N + 1:end);
b(:, 1) = flag_poss(2, 1, end - N + 1:end);
plot(t_poss, a, 'b.'), grid on;
ylim([-1.2, 2.2]);

% figure;
% missing = kt(1) * (xct_poss(475, 1) - xct_poss(474, 1)) / (xct_poss(475, 3) - xct_poss(474, 3))
% mu_min = kt(1) * (xct_poss(347, 1) - w_poss(1, 1, N + 346)) / Ft_poss(347 + N, 3)


%% wrong Numerical
clear
% A_p = load('data/Analytical Petrov System 3/kt = 30, mu = 8, k = 40, f = 100sin(t)/para_H5_ds0.05_N512_stopped_point.mat');
% A_p = load('data/Analytical Petrov System 3/kt = 30, mu = 8, k = 40, f = 100sin(t)/para_H5_ds0.05_N128_stopped_point.mat');
A_p = load('data/Analytical Petrov System 3/kt = 30, mu = 8, k = 40, f = 100sin(t)/para_H5_ds0.09_N128_stopped_point_backward_mu10.mat');

% i_plot = 707; % point before second 11.6743
x_poss = A_p.para.X(:);
% omega_poss = omega_cont(i_plot);
% omega_poss = 11.6743;
omega_poss = A_p.para.omega;
params = A_p.para.params;
xp = A_p.para.xp;
Na = A_p.para.params.func.HBM.Na;
Nx = A_p.para.params.func.HBM.Nx;
N = A_p.para.params.func.HBM.N;
H = A_p.para.params.func.HBM.H;
kn = A_p.para.params.func.fc.kn;
xn0 = A_p.para.params.func.fc.xn0;
mu = A_p.para.params.func.fc.mu;
kt = A_p.para.params.func.fc.kt;
nloop = A_p.para.params.func.fc.nloop;
for i = 1:Na + 3 * Nx
    r1 = (2 * H + 1) * (i - 1) + 1;
    r2 = (2 * H + 1) * i;
    X(:,i) = x_poss(r1:r2); % reorder in dofs in column
end
xt_poss = params.func.HBM.E * X;
xct_poss = xt_poss(:, Na + 1:end);


[FUN_poss, w_poss] = HBMFUNC(x_poss, omega_poss, params.func);
% if FUN(x) ~= 0, calculate corresponse value
if norm(FUN_poss) > params.Newton.epsf
    warning('norm(FUN_poss) > params.Newton.epsf');
    params.cont.ds = 0;
    params.cont.omega_0 = omega_poss;
    params.cont.step = 100001;
    params.cont.x0 = x_poss;
    [x_poss, omega_poss, ~, ~, w_poss, ~] = cont_step(@HBMFUNC, @HBMJACOB, @HBMJOmega, params);
    params.func.fc.w = w_poss;

    for i = 1:Na + 3 * Nx
        r1 = (2 * H + 1) * (i - 1) + 1;
        r2 = (2 * H + 1) * i;
        X(:,i) = x_poss(r1:r2); % reorder in dofs in column
    end
    xt_poss = params.func.HBM.E * X; 
    xct_poss = xt_poss(:, Na + 1:end);
end


[Ft_poss, w_poss] = g(xct_poss + xp', params.func.fc);

JNL_poss = finite_diff_jac(@(x) fftgx(x, params.func), x_poss((2 * H + 1) * Na + 1:end));

[JNL_time, Fi, xi] = finite_diff_jac_F(@(x) gFx(x, params.func), x_poss((2 * H + 1) * Na + 1:end));

T = 2 * pi / omega_poss;
dt = T / N;
t_poss = [0:dt:(T - dt)]';

para.t = t_poss;
para.xt = xt_poss;
para.Ft = Ft_poss;
para.omega = omega_poss;
para.xp = xp;

para.params = params;
para.X = X;

% solution


para.JNL_poss = JNL_poss;


figure; % displacement
subplot(2,2,1)
yyaxis left
plot(t_poss, xt_poss(:, 2) + xp(1), 'b-', 'LineWidth', 2, 'DisplayName', 'xt'), grid on;
% figure; % friction
yyaxis right
plot(t_poss, Ft_poss(end - N + 1:end, 1), 'r-', 'LineWidth', 2, 'DisplayName', 'T'), grid on;
legend('show');
titlename = 'Omega = ' + string(omega_poss) + ', Amplitude = ' + string(max(abs(xct_poss(:, 1))));
title(titlename);

subplot(2,2,2)
yyaxis left
plot(t_poss, xt_poss(:, 4) + xp(3), 'b-', 'LineWidth', 2, 'DisplayName', 'xn'), grid on;
% figure; % friction
yyaxis right
plot(t_poss, Ft_poss(end - N + 1:end, 3), 'r-', 'LineWidth', 2, 'DisplayName', 'Fn'), grid on;
legend('show');
titlename = 'Omega = ' + string(omega_poss) + ', Amplitude = ' + string(max(abs(xct_poss(:, 3))));
title(titlename);

subplot(2,2,3)
plot(xt_poss(:, 2) + xp(1), Ft_poss(end - N + 1:end, 1), 'b-', 'LineWidth', 2), grid on;
% legend('show');
xlabel('xt');
ylabel('T');

subplot(2,2,4)
plot(t_poss, mu(1) * Ft_poss(end - N + 1:end, 3), 'k-', 'LineWidth', 2, 'DisplayName', 'mu*Fn'), hold on;
plot(t_poss, -mu(1) * Ft_poss(end - N + 1:end, 3), 'k-', 'LineWidth', 2, 'DisplayName', '-mu*Fn'), grid on;
plot(t_poss, Ft_poss(end - N + 1:end, 1), 'b-', 'LineWidth', 2, 'DisplayName', 'T');
legend('show');
xlabel('t');
ylabel('F');
% savename = 'data/Analytical Petrov System 1/ky = 120, g = 10/Omega = ' + string(omega_poss) + ', Amplitude = ' + string(Adof(i_plot, 5)) + '.mat';
% save(savename, 'para');

% % plot Fi⁺ and JNL in time
% for i = 1:11 % JNL in time
%     figure;
%     if i == 1
%         titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'cos0';
%     elseif mod(i, 2) == 0
%         titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'cos' + string(floor(i/2));
%     else
%         titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'sin' + string(floor(i/2));
%     end
%     subplot(2,2,1)
%     plot(t_poss, JNL_time(1:128, i), 'k:', 'LineWidth', 2), grid on, hold on;
%     title(titlename);
%     subplot(2,2,2)
%     plot(t_poss, JNL_time(1:128, i + 22), 'k:', 'LineWidth', 2), grid on, hold on;
%     title(titlename);
%     subplot(2,2,3)
%     plot(t_poss, JNL_time(257:end, i), 'k:', 'LineWidth', 2), grid on, hold on;
%     title(titlename);
%     subplot(2,2,4)
%     plot(t_poss, JNL_time(257:end, i + 22), 'k:', 'LineWidth', 2), grid on, hold on;
%     title(titlename);
% 
% end


% % Fi
% figure;
% plot(t_poss, Fi(1:128, 1), '.'), hold on, grid on;
% for i = 1:11
%     plot(t_poss, Fi(1:128, i + 1), '.'), hold on, grid on;
% end
% figure;
% plot(t_poss, Fi(1:128, 1), '.'), hold on, grid on;
% for i = 1:11
%     plot(t_poss, Fi(1:128, i + 23), '.'), hold on, grid on;
% end
% figure;
% plot(t_poss, Fi(257:end, 1), '.'), hold on, grid on;
% for i = 1:11
%     plot(t_poss, Fi(257:end, i + 23), '.'), hold on, grid on;
% end

figure;
subplot(2,2,1)
plot(t_poss, Fi(1:128, 1), 'b.', 'DisplayName','T1 original'), hold on, grid on;
plot(t_poss, Fi(1:128, 25), 'r.', 'DisplayName','T1 dcos1', 'MarkerSize', 3), hold on, grid on;
legend('show');
subplot(2,2,2)
plot(t_poss, Fi(1:128, 25) - Fi(1:128, 1), 'b.', 'DisplayName','diff'), hold on, grid on;
legend('show');
subplot(2,2,3)
plot(t_poss, xi(1:128, 1), 'b.', 'DisplayName','xt original'), hold on, grid on;
plot(t_poss, xi(1:128, 25), 'r.', 'DisplayName','xt dcos1', 'MarkerSize', 3), hold on, grid on;
legend('show');
subplot(2,2,4)
plot(t_poss, xi(1:128, 25) - xi(1:128, 1), 'b.', 'DisplayName','diff'), hold on, grid on;
legend('show');

figure;
subplot(2,2,1)
plot(t_poss, Fi(257:end, 1), 'b.', 'DisplayName','Fn original'), hold on, grid on;
plot(t_poss, Fi(257:end, 25), 'r.', 'DisplayName','Fn dcos1', 'MarkerSize', 3), hold on, grid on;
legend('show');
subplot(2,2,2)
plot(t_poss, Fi(257:end, 25) - Fi(257:end, 1), 'b.', 'DisplayName','diff'), hold on, grid on;
legend('show');
subplot(2,2,3)
plot(t_poss, xi(257:end, 1), 'b.', 'DisplayName','xn original'), hold on, grid on;
plot(t_poss, xi(257:end, 25), 'r.', 'DisplayName','xn dcos1', 'MarkerSize', 3), hold on, grid on;
legend('show');
subplot(2,2,4)
plot(t_poss, xi(257:end, 25) - xi(257:end, 1), 'b.', 'DisplayName','diff'), hold on, grid on;
legend('show');


figure;
titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'cos1';
subplot(2,2,1)
plot(t_poss, JNL_time(1:128, 2), 'k:', 'LineWidth', 2), grid on, hold on;
title(titlename);
subplot(2,2,2)
plot(t_poss, JNL_time(1:128, 2 + 22), 'k:', 'LineWidth', 2), grid on, hold on;
title(titlename);
subplot(2,2,3)
plot(t_poss, JNL_time(257:end, 2), 'k:', 'LineWidth', 2), grid on, hold on;
title(titlename);
subplot(2,2,4)
plot(t_poss, JNL_time(257:end, 2 + 22), 'k:', 'LineWidth', 2), grid on, hold on;
title(titlename);
%% Jacobian Comparison
clear
para_A = load('data/Analytical Petrov System 3/kt = 30, mu = 8, k = 40, f = 100sin(t)/para_H5_ds0.09_N128_stopped_point_backward_mu10.mat');
para_N = load('data/Numerical Petrov System 3/kt = 30, mu = 8, k = 40, f = 100sin(t)/para_H5_ds0.09_N128_stopped_point_backward_mu10.mat');
para_F = load('data/Fixed Numerical Petrov System 3/kt = 30, mu = 8, k = 40, f = 100sin(t)/para_H5_ds0.09_N128_stopped_point_backward_mu10.mat');
E = para_A.para.params.func.HBM.E;
H = para_A.para.params.func.HBM.H;
N = para_A.para.params.func.HBM.N;

figure
yname = 'CB1';
for i = 1:4
    subplot(2, 2, i)
    plot(para_N.para.t, para_N.para.xt(:, i), 'b-', 'LineWidth', 2), hold on;
    plot(para_A.para.t, para_A.para.xt(:, i), 'k:', 'LineWidth', 2), grid on;
    plot(para_F.para.t, para_F.para.xt(:, i), 'r--', 'LineWidth', 2), grid on;
    if i > 1
        yname = 'x' + string(i - 1);
    end
    ylabel(yname);
    legend('Numerical', 'Analytical', 'Fixed Numerical');
end

figure; % displacement & forces & cycle
subplot(2,3,1)
plot(para_N.para.t, para_N.para.xt(:, 2), 'k:', 'LineWidth', 2, 'DisplayName', 'Numerical'), hold on;
plot(para_A.para.t, para_A.para.xt(:, 2), 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical'), grid on;
plot(para_F.para.t, para_F.para.xt(:, 2), 'r--', 'LineWidth', 2, 'DisplayName', 'Fixed Numerical'), legend('show');
title('Displacement x1');
subplot(2,3,2)
plot(para_N.para.t, para_N.para.xt(:, 4), 'k:', 'LineWidth', 2, 'DisplayName', 'Numerical'), hold on;
plot(para_A.para.t, para_A.para.xt(:, 4), 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical'), grid on;
plot(para_F.para.t, para_F.para.xt(:, 4), 'r--', 'LineWidth', 2, 'DisplayName', 'Fixed Numerical'), legend('show');
title('Displacement xn');
subplot(2,3,3)
plot(para_N.para.xt(:, 2), para_N.para.Ft(end - N + 1:end, 1), 'k:', 'LineWidth', 2, 'DisplayName', 'Numerical'), hold on;
plot(para_A.para.xt(:, 2), para_A.para.Ft(end - N + 1:end, 1), 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical'), grid on;
plot(para_F.para.xt(:, 2), para_F.para.Ft(end - N + 1:end, 1), 'r--', 'LineWidth', 2, 'DisplayName', 'Fixed Numerical'), legend('show');
title('Ft VS x1');
subplot(2,3,4)
plot(para_N.para.t, para_N.para.Ft(end - N + 1:end, 1), 'k:', 'LineWidth', 2, 'DisplayName', 'Numerical'), hold on;
plot(para_A.para.t, para_A.para.Ft(end - N + 1:end, 1), 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical'), grid on;
plot(para_F.para.t, para_F.para.Ft(end - N + 1:end, 1), 'r--', 'LineWidth', 2, 'DisplayName', 'Fixed Numerical'), legend('show');
title('Ft');
subplot(2,3,5)
plot(para_N.para.t, para_N.para.Ft(end - N + 1:end, 1), 'k*', 'LineWidth', 2, 'DisplayName', 'Numerical'), hold on;
plot(para_A.para.t, para_A.para.Ft(end - N + 1:end, 1), 'bo', 'LineWidth', 2, 'DisplayName', 'Analytical'), grid on;
plot(para_F.para.t, para_F.para.Ft(end - N + 1:end, 1), 'r+', 'LineWidth', 2, 'DisplayName', 'Fixed Numerical'), legend('show');
title('Ft');
subplot(2,3,6)
plot(para_N.para.t, para_N.para.Ft(end - N + 1:end, 3), 'k:', 'LineWidth', 2, 'DisplayName', 'Numerical'), hold on;
plot(para_A.para.t, para_A.para.Ft(end - N + 1:end, 3), 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical'), grid on;
plot(para_F.para.t, para_F.para.Ft(end - N + 1:end, 3), 'r--', 'LineWidth', 2, 'DisplayName', 'Fixed Numerical'), legend('show');
title('Fn');

eps_J1 = norm(para_A.para.JNL_poss - para_N.para.JNL_poss) / norm(para_A.para.JNL_poss);
eps_J2 = norm(para_A.para.JNL_poss - para_F.para.JNL_poss) / norm(para_A.para.JNL_poss);
JNLt_A_11 = E * para_A.para.JNL_poss(1:2 * H + 1, 1:2 * H + 1);
JNLt_A_13 = E * para_A.para.JNL_poss(1:2 * H + 1, 2 * (2 * H + 1) + 1:end);
JNLt_A_31 = E * para_A.para.JNL_poss(2 * (2 * H + 1) + 1:end, 1:2 * H + 1);
JNLt_A_22 = E * para_A.para.JNL_poss((2 * H + 1) + 1:2 * (2 * H + 1), (2 * H + 1) + 1:2 * (2 * H + 1));
JNLt_A_23 = E * para_A.para.JNL_poss((2 * H + 1) + 1:2 * (2 * H + 1), 2 * (2 * H + 1) + 1:end);
JNLt_A_32 = E * para_A.para.JNL_poss(2 * (2 * H + 1) + 1:end, (2 * H + 1) + 1:2 * (2 * H + 1));
JNLt_A_33 = E * para_A.para.JNL_poss(2 * (2 * H + 1) + 1:end, 2 * (2 * H + 1) + 1:end);

JNLt_N_11 = E * para_N.para.JNL_poss(1:2 * H + 1, 1:2 * H + 1);
JNLt_N_13 = E * para_N.para.JNL_poss(1:2 * H + 1, 2 * (2 * H + 1) + 1:end);
JNLt_N_31 = E * para_N.para.JNL_poss(2 * (2 * H + 1) + 1:end, 1:2 * H + 1);
JNLt_N_22 = E * para_N.para.JNL_poss((2 * H + 1) + 1:2 * (2 * H + 1), (2 * H + 1) + 1:2 * (2 * H + 1));
JNLt_N_23 = E * para_N.para.JNL_poss((2 * H + 1) + 1:2 * (2 * H + 1), 2 * (2 * H + 1) + 1:end);
JNLt_N_32 = E * para_N.para.JNL_poss(2 * (2 * H + 1) + 1:end, (2 * H + 1) + 1:2 * (2 * H + 1));
JNLt_N_33 = E * para_N.para.JNL_poss(2 * (2 * H + 1) + 1:end, 2 * (2 * H + 1) + 1:end);

JNLt_F_11 = E * para_F.para.JNL_poss(1:2 * H + 1, 1:2 * H + 1);
JNLt_F_13 = E * para_F.para.JNL_poss(1:2 * H + 1, 2 * (2 * H + 1) + 1:end);
JNLt_F_31 = E * para_F.para.JNL_poss(2 * (2 * H + 1) + 1:end, 1:2 * H + 1);
JNLt_F_22 = E * para_F.para.JNL_poss((2 * H + 1) + 1:2 * (2 * H + 1), (2 * H + 1) + 1:2 * (2 * H + 1));
JNLt_F_23 = E * para_F.para.JNL_poss((2 * H + 1) + 1:2 * (2 * H + 1), 2 * (2 * H + 1) + 1:end);
JNLt_F_32 = E * para_F.para.JNL_poss(2 * (2 * H + 1) + 1:end, (2 * H + 1) + 1:2 * (2 * H + 1));
JNLt_F_33 = E * para_F.para.JNL_poss(2 * (2 * H + 1) + 1:end, 2 * (2 * H + 1) + 1:end);

t = para_A.para.t;
for i = 1:5
    fig = figure; %('PaperOrientation','landscape','PaperUnits','centimeters','PaperPosition', 100 * [0 0 29.7 21], 'PaperSize',[29.7 21] * 100);
    
    fig.WindowState = 'maximized';
    if i == 1
        titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'cos0';
    elseif mod(i, 2) == 0
        titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'cos' + string(floor(i/2));
    else
        titlename = 'N = ' + string(N) + ', H = ' + string(H) + ', ' +'sin' + string(floor(i/2));
    end
    title(titlename);

    subplot(2,2,1)
    plot(t, JNLt_A_11(:, i), 'b-', 'LineWidth', 2), hold on;
    plot(t, JNLt_N_11(:, i), 'k:', 'LineWidth', 2), grid on;
    plot(t, JNLt_F_11(:, i), 'r--', 'LineWidth', 2), grid on;
    legend('dT1dxt A', 'dT1dxt N', 'dT1dxt F');
    title(titlename);
    % subplot(2,4,5)
    % % plot(t, (JNLt_A_11(:, i) - JNLt_N_11(:, i)) ./ norm(JNLt_N_11(:, i)), 'LineWidth', 2), hold on
    % plot(t, (JNLt_A_11(:, i) - JNLt_N_11(:, i)) , 'LineWidth', 2), hold on
    % % legend('relative difference');
    % legend('difference');
    % grid on

    subplot(2,2,2)
    plot(t, JNLt_A_13(:, i), 'b-', 'LineWidth', 2), hold on;
    plot(t, JNLt_N_13(:, i), 'k:', 'LineWidth', 2), grid on;
    plot(t, JNLt_F_13(:, i), 'r--', 'LineWidth', 2), grid on;
    legend('dT1dxn A', 'dT1dxn N', 'dT1dxn F');
    title(titlename);
    % subplot(2,4,6)
    % % plot(t, (JNLt_A_13(:, i) - JNLt_N_13(:, i)) ./ norm(JNLt_N_13(:, i)), 'LineWidth', 2), hold on
    % plot(t, (JNLt_A_13(:, i) - JNLt_N_13(:, i)), 'LineWidth', 2), hold on
    % % legend('relative difference');
    % legend('difference');
    % grid on

    subplot(2,2,3)
    plot(t, JNLt_A_31(:, i), 'b-', 'LineWidth', 2), hold on;
    plot(t, JNLt_N_31(:, i), 'k:', 'LineWidth', 2), grid on;
    plot(t, JNLt_F_31(:, i), 'r--', 'LineWidth', 2), grid on;
    legend('dFndxt A', 'dFndxt N', 'dFndxt F');
    title(titlename);
    % subplot(2,4,7)
    % % plot(t, (JNLt_A_22(:, i) - JNLt_N_22(:, i)) ./ norm(JNLt_N_22(:, i)), 'LineWidth', 2), hold on
    % plot(t, (JNLt_A_22(:, i) - JNLt_N_22(:, i)), 'LineWidth', 2), hold on
    % % legend('relative difference');
    % legend('difference');
    % grid on

    subplot(2,2,4)
    plot(t, JNLt_A_33(:, i), 'b-', 'LineWidth', 2), hold on;
    plot(t, JNLt_N_33(:, i), 'k:', 'LineWidth', 2), grid on;
    plot(t, JNLt_F_33(:, i), 'r--', 'LineWidth', 2), grid on;
    legend('dFndxn A', 'dFndxn N', 'dFndxn F');
    title(titlename);
    % subplot(2,4,8)
    % % plot(t, (JNLt_A_23(:, i) - JNLt_N_23(:, i)) ./ norm(JNLt_N_23(:, i)), 'LineWidth', 2), hold on
    % plot(t, (JNLt_A_23(:, i) - JNLt_N_23(:, i)), 'LineWidth', 2), hold on
    % % legend('relative difference');
    % legend('difference');
    % grid on

    % pintname = string(titlename) + '.png';
    % set(gcf,'PaperOrientation','landscape');
    % exportgraphics(gcf, pintname,'ContentType','vector');

end
