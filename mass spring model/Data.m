%% FEM matrices for test

FEM.Mee = eye(3);
FEM.Mcc = eye(3);
FEM.Mec = zeros(3,3);


FEM.Kee = diag([2, 2, 2]);
% FEM.Kee(1:3, 4:6) = diag([-1, -2, -3]);
% FEM.Kee(4:6, 1:3) = diag([-1, -2, -3]);

FEM.Kec = diag([-1, -1, -1]);

FEM.Kcc = diag([1, 1, 1 ]);

FEM.Fe = 1 * [1, 0.5, 0.1]';
FEM.Fc = [0, 0 ,0]';

FEM.Pe = [0, 0, 1]';
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

epsx = 1e-3;
epsf = 1e-3;
maxiter = 100;

%% Coulomb friction

kn = 1; % .* 1e5;
xn0 = 0;
mu = 0.5 * [0.1; 0.1]; % .* 1e5;
kt = 0.01 * [1; 1]; % .* 1e5;
nloop = 2;

%% preload initial condition for Newton method

xp0 = 0; % if no value defined here, the default value inside is 0


%% continuation parameters


ds = 0.01;
maxstep = 50000;
% omega_0 = 0.81; % 
% omega_end = 0.87;






% x0 = 0;

% tx0 = 0; % always 0 when calculating the first point so default inside the code
% tomega0 = 1; % is defined by lambda0 and lambda_end

