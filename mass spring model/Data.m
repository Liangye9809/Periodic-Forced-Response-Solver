%% FEM matrices for test

FEM.Mee = eye(3);
FEM.Mcc = eye(3);
FEM.Mec = zeros(3,3);


% FEM.Kee = 10 .* diag([2, 2, 0.7]);
% FEM.Kec = 10 .* diag([-1, -1, -0.35]);
% FEM.Kcc = -FEM.Kec;
% 
% FEM.Fe = 2 * [1, 0, 2]';
% FEM.Fc = [0, 0 ,0]';
% 
% FEM.Pe = 20 .* [0, 0, 1]';
% FEM.Pc = [0, 0 ,0]';
k1 = 3.5;
k2 = 10;
k3 = 3.5;
k4 = 10;

FEM.Kee = diag([(k2 + k4), (k2 + k4), (k1 + k3)]);
FEM.Kec = diag([-k4, -k4, -k3]);
FEM.Kcc = -FEM.Kec;

FEM.Fe = 1 * [2, 1, 4]';
FEM.Fc = [0, 0 ,0]';

FEM.Pe = 20 .* [0, 0, 1]';
FEM.Pc = [0, 0 ,0]';
% xn = - 4 * sin(sin(t)) + 1; % separation to stick
% xt = 2 * exp(cos(t + 1)) - 3; % separation to stick

%% External  forces

H_F_ext = [0, 0, 1]; % fourier coefficient of f(t)


%% HBM parameters

H = 5; % number of harmonics assumption
N = 2^8; % number of time points per force cycle
Nx = 1; % number of contact points, means having 4 * 3 = 12 dofs
Na = 1; % number of CB modes
xi = 0.1; 


%% Newton Method parameters

epsx = 1e-6;
epsf = 1e-6;
maxiter = 100;

%% Coulomb friction

kn = 20; % .* 1e5;
xn0 = 0;
mu = 0.5 * [1; 1]; % .* 1e5;
kt = 1 * [1; 1]; % .* 1e5;
nloop = 2;

%% preload initial condition for Newton method

xp0 = 0; % if no value defined here, the default value inside is 0


%% continuation parameters


ds = 0.01;
maxstep = 500000;
% omega_0 = 0.81; % 
% omega_end = 0.87;






% x0 = 0;

% tx0 = 0; % always 0 when calculating the first point so default inside the code
% tomega0 = 1; % is defined by lambda0 and lambda_end

