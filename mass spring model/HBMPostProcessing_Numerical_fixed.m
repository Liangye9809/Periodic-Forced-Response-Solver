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
% OMEGA = sqrt(omega02) .* omega_cont';
OMEGA = omega_cont';
Adof = [OMEGA, Adof, k_cont'];



figure
% yyaxis right
% stem(Adof(:, 2), k_cont'), grid on
% 
% yyaxis left
plot(Adof(:,2), Adof(:,3), 'b-', 'LineWidth', 2), hold on;
grid on;
ylabel('CB1');

title('Numerical Jacobian');


figure % X1
% yyaxis right
% stem(Adof(:, 2), k_cont'), grid on
% 
% yyaxis left
plot(Adof(:, 2), Adof(:, 4), 'b-', 'LineWidth', 2), hold on;
grid on;
ylabel('X1');

figure % X2
% yyaxis right
% stem(Adof(:, 2), k_cont'), grid on
% 
% yyaxis left
plot(Adof(:, 2), Adof(:, 5), 'b-', 'LineWidth', 2), hold on;
grid on;
ylabel('X2');


figure % Xn
% yyaxis right
% stem(Adof(:, 2), k_cont'), grid on

% yyaxis left
plot(Adof(:, 2), Adof(:, 6), 'b-', 'LineWidth', 2), hold on;
grid on;
ylabel('Xn');



% save Adof_Numerical_mex_g.mat Adof
%% calculate nonlinear forces
% close all
% check FUN(x) = 0;
% para_A = load("data/Analytical Jacobian results 1.0/para_gap_to_stick_point.mat");
% para_A = load("data/Analytical Jacobian results 1.1/para_gap_to_stick_point.mat");
% para_A = load("data/Analytical Jacobian results 3.0 big diff/para_gap_to_stick_point.mat");
para_A = load("data/Analytical Jacobian results 3.1 big diff no slip/para_gap_to_stick_point.mat");
x_poss = para_A.para.X(:);
omega_poss = para_A.para.omega;
params.func.fc.w = para_A.para.params.func.fc.w;

for i = 1:Na + 3 * Nx
    r1 = (2 * H + 1) * (i - 1) + 1;
    r2 = (2 * H + 1) * i;
    X(:,i) = x_poss(r1:r2); % reorder in dofs in column
end
xt_poss = params.func.HBM.E * X;
xct_poss = xt_poss(:, Na + 1:end);

[FUN_poss, w_poss, JL_poss, flag_poss] = HBMFUNC(x_poss, xct_poss, omega_poss, params.func);

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

% [Ft_poss, w_poss] = g(xct_poss, params.func.fc);

   EH = params.func.HBM.EH;
    N = params.func.HBM.N;
  gxp = params.func.static.preload.gxp;
   xp = params.func.static.preload.xp;
   kn = params.func.fc.kn;
  xn0 = params.func.fc.xn0;
   mu = params.func.fc.mu;
   kt = params.func.fc.kt;
 w_in = params.func.fc.w;
nloop = params.func.fc.nloop;

[Fti, wi, flagi] = g(xct_poss + xp', kn, xn0, mu, kt, w_in, nloop); 
Ft_poss = Fti(end - N + 1:end, :);

JNL_poss = finite_diff_jac(@(x) fftgx(x, xct_poss, params.func), x_poss(Na * (2 * H + 1) + 1:end));
T = 2 * pi / omega_poss;
dt = T / N;
t_poss = [0:dt:(T - dt)]';

para.t = t_poss;
para.xt = xt_poss;
para.Ft = Ft_poss;
para.omega = omega_poss;
para.xp = xp;
para.gxp = gxp;
% para.Pe = FEM.Pe;
% para.Pc = FEM.Pc;
para.params = params;
para.X = X;
% solution
para.x_cont = x_cont;
para.k_cont = k_cont;
para.omega_cont = omega_cont;
para.JNL_poss = JNL_poss;


figure; % displacement
plot(t_poss, xt_poss(:, 2) + xp(1), 'b-', 'LineWidth', 2), hold on;
plot(t_poss, xt_poss(:, 3) + xp(2), 'r-', 'LineWidth', 2), hold on;
plot(t_poss, xt_poss(:, 4) + xp(3), 'k-', 'LineWidth', 2), grid on;
legend('x1', 'x2', 'xn');

figure; % cycle
plot(xt_poss(:, 2), Ft_poss(end - N + 1:end, 1), 'b-', 'LineWidth', 2), hold on;
plot(xt_poss(:, 3), Ft_poss(end - N + 1:end, 2), 'r-', 'LineWidth', 2), hold on;
grid on;
xlabel('x');
ylabel('T');
legend('x1', 'x2');
title('hysteresis cycle')

figure; % friction
plot(t_poss,  mu(1) .* kn .* max((xt_poss(:, 4) + xp(3)), 0), 'k-', 'LineWidth', 2), hold on;
plot(t_poss, -mu(2) .* kn .* max((xt_poss(:, 4) + xp(3)), 0), 'k-', 'LineWidth', 2), grid on;
plot(t_poss, Ft_poss(end - N + 1:end, 1), 'b-', 'LineWidth', 2);
plot(t_poss, Ft_poss(end - N + 1:end, 2), 'r-', 'LineWidth', 2);
legend('mu * Fn', '-mu * Fn', 'T1', 'T2');


% save para_gap_to_stick_point.mat para;
%% compare amplitude

Adof_A = load("data/Analytical Jacobian results 2.0 simple sin/Adof_Analytical_matlab_g.mat");
Adof_N = load("data/Fixed Numerical Jacobian results 2.0 simple sin/Adof_Numerical_mex_g.mat");
figure
for i = 1:4
    subplot(2,2,i)
    plot(Adof_N.Adof(:, 2), Adof_N.Adof(:, i + 2), 'b-', 'LineWidth', 2), hold on;
    plot(Adof_A.Adof(:, 2), Adof_A.Adof(:, i + 2), 'r--', 'LineWidth', 2), grid on;
    yname = 'CB1';
    if i > 1
        yname = 'x' + string(i - 1);
    end
    ylabel(yname);
    legend('Numerical', 'Analytical');
end

%% compare iteration
figure
Adof_A = load("data/Analytical Jacobian results 2.0 simple sin/Adof_Analytical_matlab_g.mat");
Adof_N = load("data/Fixed Numerical Jacobian results 2.0 simple sin/Adof_Numerical_mex_g.mat");
para_A = load("data/Analytical Jacobian results 2.0 simple sin/para_gap_to_stick_point.mat");
plot(Adof_N.Adof(:, 2), Adof_N.Adof(:, 7), 'b*'), hold on, grid on;

plot(Adof_A.Adof(:, 2), Adof_A.Adof(:, 7)', 'ro'), hold on, grid on;
legend('Numerical', 'Analytical');

%% compare displacement and Jacobian
para_A = load("data/Analytical Jacobian results 2.0 simple sin/para_gap_to_stick_point.mat");
para_N = load('data/Fixed Numerical Jacobian results 2.0 simple sin/para_gap_to_stick_point.mat');
E = para_A.para.params.func.HBM.E;
H = para_A.para.params.func.HBM.H;
N = para_A.para.params.func.HBM.N;
figure
yname = 'CB1';
for i = 1:4
    subplot(2, 2, i)
    plot(para_N.para.t, para_N.para.xt(:, i), 'b-', 'LineWidth', 2), hold on;
    plot(para_A.para.t, para_A.para.xt(:, i), 'r--', 'LineWidth', 2), grid on;
    if i > 1
        yname = 'x' + string(i - 1);
    end
    ylabel(yname);
    legend('Numerical', 'Analytical');
end

eps_J = norm(para_A.para.JNL_poss - para_N.para.JNL_poss) / norm(para_A.para.JNL_poss);
JNLt_A_11 = E * para_A.para.JNL_poss(1:2 * H + 1, 1:2 * H + 1);
JNLt_A_13 = E * para_A.para.JNL_poss(1:2 * H + 1, 2 * (2 * H + 1) + 1:end);
JNLt_A_22 = E * para_A.para.JNL_poss((2 * H + 1) + 1:2 * (2 * H + 1), (2 * H + 1) + 1:2 * (2 * H + 1));
JNLt_A_23 = E * para_A.para.JNL_poss((2 * H + 1) + 1:2 * (2 * H + 1), 2 * (2 * H + 1) + 1:end);

JNLt_N_11 = E * para_N.para.JNL_poss(1:2 * H + 1, 1:2 * H + 1);
JNLt_N_13 = E * para_N.para.JNL_poss(1:2 * H + 1, 2 * (2 * H + 1) + 1:end);
JNLt_N_22 = E * para_N.para.JNL_poss((2 * H + 1) + 1:2 * (2 * H + 1), (2 * H + 1) + 1:2 * (2 * H + 1));
JNLt_N_23 = E * para_N.para.JNL_poss((2 * H + 1) + 1:2 * (2 * H + 1), 2 * (2 * H + 1) + 1:end);

t = para_A.para.t;
for i = 1:11
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

    subplot(2,4,1)
    plot(t, JNLt_A_11(:, i), 'LineWidth', 2), hold on
    plot(t, JNLt_N_11(:, i), 'LineWidth', 2), grid on;
    legend('dTdxt A', 'dTdxt N');
    title(titlename);
    subplot(2,4,5)
    % plot(t, (JNLt_A_11(:, i) - JNLt_N_11(:, i)) ./ norm(JNLt_N_11(:, i)), 'LineWidth', 2), hold on
    plot(t, (JNLt_A_11(:, i) - JNLt_N_11(:, i)) , 'LineWidth', 2), hold on
    % legend('relative difference');
    legend('difference');
    grid on

    subplot(2,4,2)
    plot(t, JNLt_A_13(:, i), 'LineWidth', 2), hold on
    plot(t, JNLt_N_13(:, i), 'LineWidth', 2), grid on;
    legend('dTdxn A', 'dTdxn N');
    title(titlename);
    subplot(2,4,6)
    % plot(t, (JNLt_A_13(:, i) - JNLt_N_13(:, i)) ./ norm(JNLt_N_13(:, i)), 'LineWidth', 2), hold on
    plot(t, (JNLt_A_13(:, i) - JNLt_N_13(:, i)), 'LineWidth', 2), hold on
    % legend('relative difference');
    legend('difference');
    grid on

    subplot(2,4,3)
    plot(t, JNLt_A_22(:, i), 'LineWidth', 2), hold on
    plot(t, JNLt_N_22(:, i), 'LineWidth', 2), grid on;
    legend('dTdxn A', 'dTdxn N');
    title(titlename);
    subplot(2,4,7)
    % plot(t, (JNLt_A_22(:, i) - JNLt_N_22(:, i)) ./ norm(JNLt_N_22(:, i)), 'LineWidth', 2), hold on
    plot(t, (JNLt_A_22(:, i) - JNLt_N_22(:, i)), 'LineWidth', 2), hold on
    % legend('relative difference');
    legend('difference');
    grid on

    subplot(2,4,4)
    plot(t, JNLt_A_23(:, i), 'LineWidth', 2), hold on
    plot(t, JNLt_N_23(:, i), 'LineWidth', 2), grid on;
    legend('dTdxn A', 'dTdxn N');
    title(titlename);
    subplot(2,4,8)
    % plot(t, (JNLt_A_23(:, i) - JNLt_N_23(:, i)) ./ norm(JNLt_N_23(:, i)), 'LineWidth', 2), hold on
    plot(t, (JNLt_A_23(:, i) - JNLt_N_23(:, i)), 'LineWidth', 2), hold on
    % legend('relative difference');
    legend('difference');
    grid on

    % pintname = string(titlename) + '.png';
    % set(gcf,'PaperOrientation','landscape');
    % exportgraphics(gcf, pintname,'ContentType','vector');

end


