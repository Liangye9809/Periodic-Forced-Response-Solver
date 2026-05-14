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

ind_gap_stick = find(gap_cont == 1 & (slipP_cont + slipM_cont) == 0);
ind_gap = find(gap_cont == 1);
ind_slip = find(slipP_cont == 1);

figure;
% yyaxis left
plot(Adof(:, 1), Adof(:, 5), 'r-', 'LineWidth', 2), hold on;
plot(Adof(:, 1), Adof(:, 3), 'b-.', 'LineWidth', 2), grid on;
xlabel('Omega');
legendname = 'g = ' + string(xn0);
legend(legendname);
% ylim([0.1, 100]);

% yyaxis right
% plot(Adof(:, 1), gap_cont', 'LineWidth', 2, 'LineStyle', '-', 'Color', 'g'), hold on;
% plot(Adof(:, 1), slipP_cont', 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'), hold on;
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



% figure;
% hold on;
% 
% % Plot with colors, line styles, and markers
% plot(x, y1, '.-',  'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'sin(x)'); % Solid
% 
% % Add labels, legend, and grid
% xlabel('x');
% ylabel('y');
% legend('show', 'Location', 'best');
% grid on;
% hold off;

% figure;
% plot(Adof(:, 1), Adof(:, 5), 'r-', 'LineWidth', 2), hold on;
% grid on;



figure
hold on;
load('data/Analytical Petrov System 1/ky = 120, g = 10000000000.mat');
plot(Adof(:, 1), Adof(:, 5), 'LineWidth', 2, 'DisplayName', 'separated'), hold on;

load('data/Analytical Petrov System 1/ky = 120, full contact.mat');
plot(Adof(:, 1), Adof(:, 5), 'LineWidth', 2, 'DisplayName', 'full contact'), hold on;
for i = -10:5:10
    filemane = 'data/Analytical Petrov System 1/ky = 120, g = ' + string(i) + '.mat';
    load(filemane);
    plot(Adof(:, 1), Adof(:, 5), 'LineWidth', 2, 'DisplayName', string(i)), hold on;
end
xlabel('Frequency, rad/s');
ylabel('Maximum displacement');
legend('show');
grid on;

%%


i_plot = 2010;
x_poss = x_cont(:, i_plot);
omega_poss = omega_cont(i_plot);
params.func.fc.w = w_cont(:, i_plot);

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
para.x_cont = x_cont;
para.k_cont = k_cont;
para.omega_cont = omega_cont;
para.slipM_cont = slipM_cont;
para.slipP_cont = slipP_cont;
para.stick_cont = stick_cont;
para.gap_cont = gap_cont;
para.JNL_poss = JNL_poss;


figure; % displacement
yyaxis left
plot(t_poss, xt_poss(:, 4) + xp(3), 'b-', 'LineWidth', 2, 'DisplayName', 'xn'), grid on;


% figure; % friction
yyaxis right
plot(t_poss, Ft_poss(end - N + 1:end, 3), 'r-', 'LineWidth', 2, 'DisplayName', 'Fn'), grid on;
legend('show');

titlename = 'Omega = ' + string(omega_poss) + ', Amplitude = ' + string(Adof(i_plot, 5));
title(titlename);
savename = 'data/Analytical Petrov System 1/ky = 120, g = 10/Omega = ' + string(omega_poss) + ', Amplitude = ' + string(Adof(i_plot, 5)) + '.mat';
save(savename, 'para');