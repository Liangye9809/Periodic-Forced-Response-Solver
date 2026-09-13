%% mu = 0.5 case
% pathin = pwd;
% % cd('/home/liangye-liu/data/non-linear problem/Periodic-Forced-Response-Solver/mass spring model/data/Analytical J and F');
% cd('D:\study\PHD\data\Frictions\Periodic-Forced-Response-Solver\mass spring model\data\Analytical J and F');
% % load("Data_unconverge.mat");
% load("Data_unconverge_mu0.5_slip_gap_stick_1.mat");
% cd(pathin);
% X = D.x;
% xp = D.params.func.static.preload.xp;
% gxp = D.params.func.static.preload.gxp;
% % S = D.S;
% N = D.params.func.HBM.N;
% H = D.params.func.HBM.H;
% E = D.params.func.HBM.E;
% EH = D.params.func.HBM.EH;
% xct = D.xct + xp';
% kt = D.params.func.fc.kt;
% kn = D.params.func.fc.kn;
% mu = D.params.func.fc.mu;
% w = D.params.func.fc.w;
% nloop = D.params.func.fc.nloop;
% Nx = D.params.func.HBM.Nx;
% Na = D.params.func.HBM.Na;
% xn0 = D.params.func.fc.xn0;
% 
% Xc = X(Na * (2 * H + 1) + 1:end);
% for ixp = 1:3 * Nx
%     Xc((2 * H + 1) * (ixp - 1) + 1) = Xc((2 * H + 1) * (ixp - 1) + 1) + 2 * xp(ixp); % consider preload
% end


%% simple case
N = 2^12;
H = 100;
dt = 2 * pi / N;
t = (0:(N-1)) * 2 * pi / N;
t = t';
% xn = ones(N, 1);
% xt = 0.5 * sin(t); % pure stick
% xn = - 4 * sin(sin(t)) + 1; % separation to stick
% xt = 2 * exp(cos(t + 1)) - 3; % separation to stick
xn = 2 * exp(cos(t)) - 0.5; % slip to stick
% xn = 2 * exp(cos(t)) - 0.75; % separation to slip
xt = 2 * sin(sin(t)); % slip to stick

% xt = 1.05  * sin(2 .* exp(cos(t))); % tangent case1
% xt = 1.00 * sin(sin(t)) ./ sin(1); % tangent case2
% xt = 0.5 * sin(sin(t)) ./ sin(1) + 0.5; % tangent case3 only one side

% plot w for the paper (kt = 1, kn = 2, mu = 0.5)
% pure stick
% xn = 2*ones(N, 1); % pure stick
% xt = sin(sin(t)); % pure stick
% simple x
% xn = 2*ones(N, 1); % pure stick
% xt = sin(t); % pure stick

% slip to stick
% xn = 2 * exp(cos(t)) - 0.5; % slip to stick
% xt = 2 * sin(sin(t)); % slip to stick
% simple x
% xn = 2.5 * cos(t) + 3; % slip to stick
% xt = 3 * sin(t); % slip to stick

% gap to stick
% xn = - 4 * sin(sin(t)) + 1; % separation to stick
% xt = 2 * exp(cos(t + 1)) - 3; % separation to stick
% simple x
% xn = - 10 * cos(t) + 3; % separation to stick
% xt = - sin(t); % separation to stick
% xn = - 1 * cos(t) + 0.5; % separation to stick % (kt = 1, kn = 500, mu = 0.8)
% xt = 2 * cos(t); % separation to stick
xct_ = [xt, xt, xn];
xp = [0; 0; 0];
gxp= zeros(3, 1);
Nx = 1;
Na = 0;



kt = [1; 1];
kn = 2;
mu = 0.5 * [1; 1];
w =  0 * [1; 1];
xn0 = 0; % normal pre-displacement

nloop = 2;
[E, EH] = fft_matrices(N, H);
Xc = EH * xct_;
xct = E * Xc;
Xc = Xc(:);