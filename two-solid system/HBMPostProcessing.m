%% plot Amplitude vs omega
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultAxesFontSize',16)
set(0, 'DefaultFigurePosition', [500, 500, 600, 450]);
set(0, 'DefaultFigureColor', 'w');


%%
[E, EH] = HBM.fft_matrices(2^12, H);
for i = 1:3*Nx + Na
    x_contNx(:,:,i) = x_cont((2*H+1)*(i-1)+1:(2*H+1)*i,:);
end
Adof = omega_cont';
for Ndof = 1:3*Nx + Na
    A = E * x_contNx(:,:,Ndof);
    Amax = max(abs(A));
    Adof(:, Ndof + 1) = Amax';
end
OMEGA = sqrt(omega02) .* omega_cont';
Adof = [OMEGA, Adof];

figure;
plot(Adof(:,1), Adof(:, 3:Na + 2), 'LineWidth', 2), hold on, grid on;
legend('CB1','CB2','CB3','CB4','CB5','CB6','CB7','CB8','CB9');
for i = 1:Nx
    ind_gap_stick{i} = find(gap_cont(i, :) == 1 & (slipP_cont(i, :) + slipM_cont(i, :)) == 0);
    ind_gap{i} = find(gap_cont(i, :) == 1);
    ind_slip{i} = find(slipP_cont(i, :) == 1);
end
fig1 = figure;
% yyaxis left
plot(Adof(:,1), Adof(:,Na + 3:3:end), 'LineWidth', 2), hold on, grid on;
for i = 1:Nx
    plot(Adof(ind_slip{i},1), Adof(ind_slip{i}, 3 * i + Na), 'ko', 'LineWidth', 2), hold on, grid on;
end
legend('Nx1','Nx2','Nx3','Nx4','Nx1_slip','Nx2_slip','Nx3_slip','Nx4_slip');
title('$x_1$');

fig2 = figure;
% yyaxis left
plot(Adof(:,1), Adof(:,Na + 4:3:end), 'LineWidth', 2), hold on, grid on;
for i = 1:Nx
    plot(Adof(ind_slip{i},1), Adof(ind_slip{i}, 3 * i + Na + 1), 'ko', 'LineWidth', 2), hold on, grid on;
end
legend('Nx1','Nx2','Nx3','Nx4','Nx1_slip','Nx2_slip','Nx3_slip','Nx4_slip');
title('$x_2$');

fign = figure;
% yyaxis left
plot(Adof(:,1), Adof(:,Na + 5:3:end), 'LineWidth', 2), hold on, grid on;
for i = 1:Nx
    plot(Adof(ind_gap{i},1), Adof(ind_gap{i}, 3 * i + Na + 2), 'ko', 'LineWidth', 2), hold on, grid on;
end
legend('Nx1','Nx2','Nx3','Nx4','Nx1_gap','Nx2_gap','Nx3_gap','Nx4_gap');
title('$x_n$');

% yyaxis right
% stem(Adof(:, 2), k_cont'), grid on

% title('Numerical Jacobian');
% title('Analytical Jacobian');

% Ndof = 1;
% A = E * x_contNx(:,:,Ndof);
% Amax = max(abs(A));
% figure;
% frq = omega_cont; % ./ (2*pi);
% % frq = omega_cont ./ (2*pi);
% plot(frq, Amax,'r--'),hold on;
% % xlabel('Omega [Hz]');
% xlabel('Omega');
% ylabel('||dof_' + string(Ndof) + '(t)||');
% % ylabel('Max Amplitude of x_3');
% legend('H = ' + string(H));
% grid on;
% a1 = [frq', Amax'];

% frq = omega_cont; 
% a1 = frq';
% for Ndof = [1, Na + 1, Na + 2, Na + 3]
%     A = E * x_contNx(:,:,Ndof);
%     Amax = max(abs(A));
% 
%     a1 = [a1, Amax']; % frq, a1, x1, x2, x3
% end
% save('data/NewMesh/dof_H' + string(H),'a1');
% save('data/Cont_Friction_F1H_a1_H' + string(H),'a1');

%% calculate nonlinear forces

% n_plot = 1167; % gap begin
% n_plot = 1307; % gap max
n_plot = 1447; % gap end
omega_poss = omega_cont(n_plot);
x_poss = x_cont(:, n_plot);
w_poss = w_cont(:, n_plot*4 - 3:n_plot*4);
[Ft_poss, wt_poss, flag_poss, xct_poss, dxct_poss, t_poss, Xc] = PlotOneResult(x_poss, w_poss, omega_poss, params);


para.t = t_poss;
para.xct = xct_poss;
para.Ft = Ft_poss;
para.omega = omega_poss;
para.xp = xp;
para.gxp = gxp;
para.fc = params.func.fc;
para.HBM = params.func.HBM;
% para.Pe = FEM.Pe;
% para.Pc = FEM.Pc;
para.params = params;
para.X = x_poss;

figure(101);
plot(t_poss, xct_poss(:, 1:3:end) + xp(1:3:end)', 'LineWidth', 2), hold on, grid on;
legend('Nx1','Nx2','Nx3','Nx4');
title('x1');
ylabel('Displacement');
xlabel('Omega');

figure(102);
plot(t_poss, xct_poss(:, 2:3:end) + xp(2:3:end)', 'LineWidth', 2), hold on, grid on;
legend('Nx1','Nx2','Nx3','Nx4');
title('x2');
ylabel('Displacement');
xlabel('Omega');

figure(103);
plot(t_poss, xct_poss(:, 3:3:end) + xp(3:3:end)', 'LineWidth', 2), hold on, grid on;
legend('Nx1','Nx2','Nx3','Nx4');
title('xn');
ylabel('Displacement');
xlabel('Omega');

figure(100);
semilogy([0:2*H]', abs(Xc));

figure(201); % Ft1
for i = 1:Nx
    subplot(2,2,i)
    plot(t_poss, mu(1) * Ft_poss(end - N + 1:end, 3 * i), 'k-', 'LineWidth', 2), hold on, grid on;
    plot(t_poss, -mu(1) * Ft_poss(end - N + 1:end, 3 * i), 'k-', 'LineWidth', 2), hold on, grid on;
    plot(t_poss, Ft_poss(end - N + 1:end, 3 * i - 2), 'LineWidth', 2), hold on, grid on;
    titlename = 'Ft1 Nx' + string(i);
    title(titlename);
    ylabel('Displacement');
    xlabel('Omega');
end
    
figure(202); % Ft2
for i = 1:Nx
    subplot(2,2,i)
    plot(t_poss, mu(1) * Ft_poss(end - N + 1:end, 3 * i), 'k-', 'LineWidth', 2), hold on, grid on;
    plot(t_poss, -mu(1) * Ft_poss(end - N + 1:end, 3 * i), 'k-', 'LineWidth', 2), hold on, grid on;
    plot(t_poss, Ft_poss(end - N + 1:end, 3 * i - 1), 'LineWidth', 2), hold on, grid on;
    titlename = 'Ft2 Nx' + string(i);
    title(titlename);
    ylabel('Displacement');
    xlabel('Omega');
end

figure(203);
plot(t_poss, Ft_poss(end - N + 1:end, 3:3:end), 'LineWidth', 2), hold on, grid on;
legend('Nx1','Nx2','Nx3','Nx4');
title('Fn');
ylabel('Displacement');
xlabel('Omega');
% save('data/Mesh32x32/Pe100eachData_PeFixed_Omega_'+ string(omega_plot) + '_Nx_' + string(Nx) + '_H' + string(H) + '_N' + string(N) + '.mat','para');
%% H validation
% figure;
% c = jet(15);
% for i = 1:2:15
%     name = 'data/NewMesh/dof_H' + string(i) + '.mat';
%     load(name);
%     frq = a1(:, 1);
% 
%     subplot(2, 2, 1); % a1
%     plot(frq, a1(:, 2), 'Color', c(i, :)), hold on;
% 
%     subplot(2, 2, 3); % t1
%     plot(frq, a1(:, 3), 'Color', c(i, :)), hold on;
% 
%     subplot(2, 2, 4); % t2
%     plot(frq, a1(:, 4), 'Color', c(i, :)), hold on;
% 
%     subplot(2, 2, 2); % n
%     plot(frq, a1(:, 5), 'Color', c(i, :)), hold on;
% 
% end
% subplot(2, 2, 1); % a1
% title('a1');
% xlabel('Omega');
% ylabel('||dof_1||');
% grid on;
% legend('H1', 'H3', 'H5', 'H7', 'H9', 'H11', 'H13', 'H15');
% 
% subplot(2, 2, 3); % t1
% plot(frq, a1(:, 3), 'Color', c(i, :)), hold on;
% title('x1-t1');
% xlabel('Omega');
% ylabel('||dof_6||');
% grid on;
% legend('H1', 'H3', 'H5', 'H7', 'H9', 'H11', 'H13', 'H15');
% 
% subplot(2, 2, 4); % t2
% plot(frq, a1(:, 4), 'Color', c(i, :)), hold on;
% title('x1-t2');
% xlabel('Omega');
% ylabel('||dof_7||');
% grid on;
% legend('H1', 'H3', 'H5', 'H7', 'H9', 'H11', 'H13', 'H15');
% 
% subplot(2, 2, 2); % n
% plot(frq, a1(:, 5), 'Color', c(i, :)), hold on;
% title('x1-n');
% xlabel('Omega');
% ylabel('||dof_8||');
% grid on;
% legend('H1', 'H3', 'H5', 'H7', 'H9', 'H11', 'H13', 'H15');
%% ifft and convert a to Xe
% Nc = 3*Nx;
% clear x_contDOF xt A_cont A_a A_x A_Xe at Xet
% [E,EH] = HBM.fft_matrices(N, H);
% for i = 1:Na % for elastic part, convert to Xe
%     x_contDOF(:,:,i) = x_cont((2*H+1)*(i-1)+1:(2*H+1)*i, :);
%     at(:,:,i) = E * x_contDOF(:,:,i);
% end
% j = 1;
% for i = (Na + 1):(Na + Nc) % for contact part
%     x_contDOF(:,:,i) = x_cont((2*H+1)*(i-1)+1:(2*H+1)*i, :);
%     xt(:,:,j) = E * x_contDOF(:,:,i);
%     j = j + 1;
% end
% at = permute(at,[3,1,2]);
% xt = permute(xt,[3,1,2]);
% xp = params.func.static.preload.xp;
% Xct = xt + xp;
% for j = 1:size(at,3)
%     Xet(:,:,j) = params.func.CBmods.Phi * at(:,:,j) + params.func.CBmods.Psi * Xct(:,:,j) + params.func.static.preload.xe0;
% end
% Xet = permute(Xet,[2,3,1]);
% Xct = permute(Xct,[2,3,1]);
% xt = permute(xt,[2,3,1]);
% 
% for i = 1:Na % for elastic part
%     A_cont(i,:) = max( abs(Xet(:,:,i) ));
% end
% j = 1;
% for i = (Na + 1):(Na + Nc) % for contact part
%     A_cont(i,:) = max(abs(xt(:,:,j)));
%     j = j + 1;
% end
% A_cont(Na + Nc + 1,:) = omega_cont;
% A_cont = A_cont';
% 


%% plot w
figure;
for i = 1:4
    subplot(2,2,i)
    yyaxis left
    plot(omega_cont, w_cont(1, i:4:end), 'r--', 'LineWidth', 2), hold on;
    plot(omega_cont, w_cont(2, i:4:end), 'b-', 'LineWidth', 2), grid on;
    ylabel('w');

    yyaxis right
    plot(omega_cont, ss_cont(i, :), 'LineWidth', 2);
    ylim([0, 1]);
    xlabel('Omega');
    ylabel('slip range');
    legend('w1', 'w2', 'slip range');
    titlename = 'Nx = ' + string(i);
    title(titlename);
end