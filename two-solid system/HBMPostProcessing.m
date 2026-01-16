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
OMEGA = sqrt(omega02) .* omega_cont';
Adof = [OMEGA, Adof];
% figure
% plot(Adof(:,1), Adof(:,3)), hold on;
% filename = 'data/Mesh32x32/Pe100each_Adof_CP' + string(Nx) + '_PreloadFixed_9points_H' + string(H) + '_N' + string(N) + '.mat';
% filename = 'data/NewMesh/Pe112.5each_Adof_CP' + string(Nx) + '_PreloadFixed_8points.mat';
% save(filename, "Adof");

% Ndof = 1;
% A = E * x_contNx(:,:,Ndof);
% Amax = max(abs(A));
% % figure;
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

% check FUN(x) = 0;
omega_poss = omega_end;
x_poss = x_cont(:, end);
[FUN_poss, w_poss] = HBMFUNC(x_poss, omega_poss, params.func);

% if FUN(x) ~= 0, calculate corresponse value
if norm(FUN_poss) > params.Newton.epsf
    params.cont.ds = 0;
    params.cont.step = 100001;
    [x_poss, omega_poss, ~, ~, w_poss] = cont_step(@HBMFUNC, @HBMJACOB, @HBMJOmega, params);
end

for i = 1:Na + 3 * Nx
    r1 = (2 * H + 1) * (i - 1) + 1;
    r2 = (2 * H + 1) * i;
    X(:,i) = x_poss(r1:r2); % reorder in dofs in column
end
xt_poss = params.func.HBM.E * X; 
xct_poss = xt_poss(:, Na + 1:end);
[Ft_poss, w_poss] = g(xct_poss + xp', params.func.fc);
T = 2 * pi / omega_poss;
dt = T / N;
t_poss = 0:dt:(T - dt);

para.t = t_poss';
para.xt = xt_poss;
para.Ft = Ft_poss;
para.omega = omega_poss;
para.xp = xp;
para.gxp = gxp;
para.fc = params.func.fc;
para.HBM = params.func.HBM;
para.Pe = FEM.Pe;
para.Pc = FEM.Pc;

save('data/Mesh32x32/Pe100eachData_PeFixed_Omega_'+ string(omega_plot) + '_Nx_64_H3_N256.mat','para');
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
