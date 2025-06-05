%%
H = params.func.HBM.H;
[E, EH] = HBM.fft_matrices(2^12, H);
for i = 1:3*Nx + Na
    x_contNx(:,:,i) = x_cont((2*H+1)*(i-1)+1:(2*H+1)*i,:);
end
A = E * x_contNx(:,:,1);
Amax = max(abs(A));
% figure;
frq = omega_cont ./ (2*pi);
plot(frq, Amax,'b-'),hold on;
xlabel('Omega [Hz]');
ylabel('||a_1(t)||');
legend('H = ' + string(H));
grid on;
% a1 = [frq', Amax'];
% cd data/'30 May - two-solid system'/
% dataname = 'Cont_Cubic_F1H_a1_H' + string(H) + '_0.25fext_gamma=5e16.mat';
% save(dataname,'a1')
% cd ../
% cd ../

% %% ifft and convert a to Xe
% clear x_contDOF xt A_cont A_a A_x A_Xe at Xet
% [E,EH] = HBM.fft_matrices(2^12, H);
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
% %% plot Julia results
% % omegaJ = readtable('Omega.csv');
% % AmaxJ = readtable('x.csv');
% % omegaJ = readtable('Omega_preload.csv');
% % AmaxJ = readtable('x_preload.csv');
% % omegaJ = readtable('Omega_friction.csv');
% % AmaxJ = readtable('x_friction.csv');
% % omegaJ = readtable('Omega_linear.csv');
% % AmaxJ = readtable('x_linear.csv');
% omegaJ = readtable('Omega_xn0-3.5_1H.csv');
% AmaxJ = readtable('x_xn0-3.5_1H.csv');
% figure;
% omegaJ.Properties.VariableNames{1} = 'X';
% AmaxJ.Properties.VariableNames{1} = 'Y';
% plot(omegaJ.X,AmaxJ.Y,'b-'),hold on;
% % plot(omegaJ.Var1,AmaxJ.Var1,'b-'),hold on;
% % xlim([0.2,1.8]);
% xlabel('Omega');
% ylabel('Max Amplitude of x_3');
% 
% 
% 
% %% plot matlab solution
% % figure;
% plot(A_cont(:,end),A_cont(:,5),'r-.'),hold on;
% xlabel('Omega');
% ylabel('x_3 Amptitude');
% grid on;
% legend('Julia','Matlab');
% title('preload');
