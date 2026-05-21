% t = [0:N-1]' * (2 * pi / N);
% % % xn = (- 4 * sin(sin(t)) + 1) .* 0.2; % separation to stick
% % % xt = (2 * exp(cos(t + 1)) - 3) .* 1; % separation to stick
% % xn = (-sin(t) + 0.5) .* 5; % separation to stick
% % xt = (cos(t)) .* 5; % separation to stick
% % Xn = params.func.HBM.EH * xn;
% % Xt = params.func.HBM.EH * xt;
% % Xn = Xn';
% % Xt = Xt';
% % figure;
% % plot(t, xt, 'b-'), hold on;
% % plot(t, xn, 'r-'), grid on;
% % legend('xt', 'xn');
% 
% Xn = zeros(1, 2 * H + 1);
% Xn([1,2,3]) = [1.5, 0, -1];
% Xt = zeros(1, 2 * H + 1);
% Xt([1,2,3]) = [0, -cos(0), sin(0)];
% Xt2 = zeros(1, 2 * H + 1);
% 
% xt = FEM.Fe(1) .* params.func.HBM.E * Xt';
% xn = FEM.Fe(3) .* params.func.HBM.E * Xn';
% % figure
% % plot(t, xt, 'b-', 'LineWidth', 2), hold on;
% % plot(t, xn, 'k-', 'LineWidth', 2), grid on;
% % legend('ft', 'fn');
% % for Na = 1
% Fe = diag(FEM.Fe) * [Xt; Xt2; Xn];
% Fa = CB.CBmods.Phi' * Fe / (alpha * omega02);
% Fx = CB.CBmods.Psi' * Fe .* beta / (alpha^2 * omega02);
% % Fa = CB.CBmods.Phi' * Fe;
% % Fx = CB.CBmods.Psi' * Fe;
% Fx = Fx';
% params.func.HBM.fftfa = Fa(:);
% params.func.HBM.fftfx = Fx(:);

% params.func.HBM.fftfx([2, 3]) = 100 * [sin(0.5*pi); cos(0.5*pi)];
% params.func.HBM.fftfx([2, 3]) = 500 * [cos(0.1*pi); -sin(0.1*pi)];
% params.func.HBM.fftfx([2, 3]) = 1000 * [cos(0.1*pi); -sin(0.1*pi)];
% params.func.HBM.fftfx([2, 3]) = 1000 * [0; 1];
