set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultAxesFontSize',16)
set(0, 'DefaultFigurePosition', [500, 500, 600, 450]);
set(0, 'DefaultFigureColor', 'w');

%% plot Amplitude vs omega

Np = 2^13;
[E, EH] = HBM.fft_matrices(Np, H);
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
plot(Adof(:, 1), Adof(:, 5), 'r-', 'LineWidth', 2), hold on;
plot(Adof(:, 1), Adof(:, 3), 'b-', 'LineWidth', 2), grid on;
xlabel('Omega');
titlename = 'Numerical N' + string(N) + ', H' + string(H) + ', $\epsilon$ ' + string(epsf) + ...
    ', $F=' + string(CB.CB_F.Fx(end)) + '\sin(t)$, ds' + string(ds);
title(titlename);
legendname = 'g = ' + string(Rx(end) / (40 + kn));
legend(legendname);

% figure;
% % yyaxis left
% subplot(1, 2, 1)
% plot(Adof(:, 1), Adof(:, 5), 'b-', 'LineWidth', 2), hold on;
% % plot(Adof([1042,1063],1),Adof([1042,1063],5), 'ro'), grid on;
% xlabel('Omega');
% % legend('$x_{n}$ Analytical N64','pick point')
% titlename = 'N' + string(64) + ', H' + string(H) + ', $\epsilon$ = 1e-6, $F=' + string(CB.CB_F.Fx(end)) + '\sin(t)$';
% title(titlename)
% grid on;
% % ylim([0.1, 100]);
% subplot(1, 2, 2)
% plot(Adof(:, 1), Adof(:, 5), 'b-', 'LineWidth', 2), hold on;
% % plot(Adof([1042,1063],1),Adof([1042,1063],5), 'ro'), grid on;
% xlabel('Omega');
% % legend('$x_{n}$ Analytical N64','pick point')
% titlename = 'N' + string(64) + ', H' + string(H) + ', $\epsilon$ = 1e-6, $F=' + string(CB.CB_F.Fx(end)) + '\sin(t)$';
% title(titlename)% ylim([0.1, 100]);
% grid on

% figure;
% for i = 1:3
%     plot(omega_cont, x_cont())

figure % Fourier Coefficients
plot(omega_cont, x_cont((2*H+1)*(4-1)+1:(2*H+1)*(4-1)+7,:), 'LineWidth', 2);
legend('$X_{c}^{0}$','$X_{c}^{1}$','$X_{s}^{1}$',...
    '$X_{c}^{2}$','$X_{s}^{2}$','$X_{c}^{3}$','$X_{s}^{3}$')
titlename = 'Numerical N' + string(N) + ', H' + string(H) + ', $\epsilon$ ' + string(epsf) + ...
    ', $F=' + string(CB.CB_F.Fx(end)) + '\sin(t)$, ds' + string(ds) + ' xn';
title(titlename);

figure
plot(omega_cont, x_cont((2*H+1)*(2-1)+1:(2*H+1)*(2-1)+7,:), 'LineWidth', 2);
legend('$X_{c}^{0}$','$X_{c}^{1}$','$X_{s}^{1}$',...
    '$X_{c}^{2}$','$X_{s}^{2}$','$X_{c}^{3}$','$X_{s}^{3}$')
titlename = 'Numerical N' + string(N) + ', H' + string(H) + ', $\epsilon$ ' + string(epsf) + ...
    ', $F=' + string(CB.CB_F.Fx(end)) + '\sin(t)$, ds' + string(ds) + ' xt';
title(titlename);
% figure
% semilogy(abs(x_cont((2*H+1)*(4-1)+1:(2*H+1)*(5-1),1042)), 'LineWidth', 2), hold on
% semilogy(abs(x_cont((2*H+1)*(4-1)+1:(2*H+1)*(5-1),1063)), 'LineWidth', 2), grid on
% title('Fourier modes of 2 cases with same $\Omega$')
% 
% figure
% stem(x_cont((2*H+1)*(4-1)+1:(2*H+1)*(5-1),1042), 'LineWidth', 2), hold on
% stem(x_cont((2*H+1)*(4-1)+1:(2*H+1)*(5-1),1063), 'LineWidth', 2), grid on
% title('Fourier modes of 2 cases with same $\Omega$')

P.x_cont = x_cont;
P.Adof = Adof;
P.N = N;
P.H = H;
P.k_cont = k_cont;
P.omega_cont = omega_cont;
P.epsx = epsx;
P.epsf = epsf;
P.ds = ds;
P.Rx = Rx;
%% xn with different pre gap

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

%% plot points in step region (Analytical Jacobian)
% indx = [1042,1047,1055,1062,1068]; % for H5, N64, ds = 0.05, eps 1e-6
N = 64;
% indx = [1042, 1055, 1068];% for H5, N64, ds = 0.05, eps 1e-6
% indx = [409, 410, 411];% stright H5, N64, ds = 0.05, eps 1e-6
% indx = [1041,1054,1067]; % for H5, N64, ds = 0.05, eps 1e-10
% indx = [1042,1056,1069]; % for H10, N64, ds = 0.05, eps 1e-6
indx = [1042,1063]; % for H10, N64, ds = 0.05, eps 1e-6
figure
for j = 1:size(indx, 2)
    i_plot = indx(j); % point before second 11.6743
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
    t = [0:dt:(T - dt)]';

    subplot(2,size(indx,2),j) % displacement
    plot(t, xct_poss(:, 3), 'b-', 'LineWidth', 2, 'DisplayName', '$x_n$'), grid on;
    legend show
    titlename = 'Omega = ' + string(omega_poss) + ', $max(|x_n|)=$ ' + string(max(abs(xct_poss(:, 3))));
    title(titlename);
    subplot(2,size(indx,2),j+size(indx,2)) % Force
    plot(t, Ft_poss(end - N + 1:end, 3), 'b-', 'LineWidth', 2, 'DisplayName', '$F_n$'), grid on;
    legend show
end
    

%% plot points in step region (Analytical Jacobian)
% indx = [1042,1047,1055,1062,1068]; % for H5, N64, ds = 0.05, eps 1e-6
N = 64;
% indx = [1042, 1055, 1068];% for H5, N64, ds = 0.05, eps 1e-6
% indx = [409, 410, 411];% stright H5, N64, ds = 0.05, eps 1e-6
% indx = [1041,1054,1067]; % for H5, N64, ds = 0.05, eps 1e-10
% indx = [1042,1056,1069]; % for H10, N64, ds = 0.05, eps 1e-6
indx = [1042,1063]; % for H10, N64, ds = 0.05, eps 1e-6
figure
for j = 1:size(indx, 2)
    i_plot = indx(j); % point before second 11.6743
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
    t = [0:dt:(T - dt)]';

    subplot(2,3,j) % displacement
    plot(t, xct_poss(:, 3), 'b-', 'LineWidth', 2, 'DisplayName', '$x_n$'), grid on;
    legend show
    titlename = 'Omega = ' + string(omega_poss) + ', $max(|x_n|)=$ ' + string(max(abs(xct_poss(:, 3))));
    title(titlename);

    subplot(2,3,3) % displacement
    plot(t, xct_poss(:, 3), 'b-', 'LineWidth', 2, 'DisplayName', '$x_n$'), grid on; hold on;
    legend show
    titlename = 'Omega = ' + string(omega_poss) + ', $max(|x_n|)=$ ' + string(max(abs(xct_poss(:, 3))));
    title(titlename);

    subplot(2,3,j+3) % Force
    plot(t, Ft_poss(end - N + 1:end, 3), 'b-', 'LineWidth', 2, 'DisplayName', '$F_n$'), grid on;
    legend show

    subplot(2,3,6) % Force
    plot(t, Ft_poss(end - N + 1:end, 3), 'b-', 'LineWidth', 2, 'DisplayName', '$F_n$'), grid on; hold on;
    legend show
end
    