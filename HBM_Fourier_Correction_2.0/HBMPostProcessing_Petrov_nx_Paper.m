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

ind_gap_stick = find(gap_cont == 1 & (slipP_cont + slipM_cont) == 0);
ind_gap = find(gap_cont == 1);
ind_slip = find(slipP_cont == 1);
figure;
% yyaxis left
plot(Adof(:, 1), Adof(:, 5), 'r-', 'LineWidth', 2, 'DisplayName', 'kn'), hold on;
plot(Adof(:, 1), Adof(:, 3), 'b-', 'LineWidth', 2, 'DisplayName', 'kt'), grid on;
plot(Adof(ind_slip, 1), Adof(ind_slip, 3), 'ko', 'DisplayName', 'slip points');
% yyaxis right
% stem(Adof(:, 2), k_cont'), grid on
% plot(Adof(:, 1), gap_cont', 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r'), hold on;
% plot(Adof(:, 1), slipP_cont', 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k'), hold on;
% ylim([0, 1.1]);

% legend('CB1', 'gap appears area', 'slip appears area');
xlabel('Omega');
titlename = 'Numerical N' + string(N) + ', H' + string(H) + ', $\epsilon$ ' + string(epsf) + ...
    ', $F=' + string(CB.CB_F.Fx(end)) + '\sin(t)$, ds' + string(ds);
title(titlename);
% legendname = 'g = ' + string(Rx(end) / (40 + kn));
% legend(legendname);

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
plot(omega_cont, x_cont((2*H+1)*(4-1)+1:(2*H+1)*(4-1)+3,:), 'LineWidth', 2);
legend('$X_{c}^{0}$','$X_{c}^{1}$','$X_{s}^{1}$',...
    '$X_{c}^{2}$','$X_{s}^{2}$','$X_{c}^{3}$','$X_{s}^{3}$')
titlename = 'Numerical N' + string(N) + ', H' + string(H) + ', $\epsilon$ ' + string(epsf) + ...
    ', $F=' + string(CB.CB_F.Fx(end)) + '\sin(t)$, ds' + string(ds) + ' xn';
title(titlename);

figure
plot(omega_cont, x_cont((2*H+1)*(2-1)+1:(2*H+1)*(2-1)+3,:), 'LineWidth', 2);
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
P.slipM_cont = slipM_cont;
P.slipP_cont = slipP_cont;
P.gap_cont = gap_cont;
P.stick_cont = stick_cont;

