
%%
% tep = pwd;
% cd FEM/
% load('DATA.mat');
% FEM = DATA;
% clear DATA
% cd(tep)

% FEM.Pe = - FEM.Pe;
%% External  forces

H_F_ext = [0, 0, 1]; % fourier coefficient of f(t)


%% HBM parameters

% H = 7; % number of harmonics assumption
% N = 2^6; % number of time points per force cycle
% Nx = 64; % number of contact points, means having 4 * 3 = 12 dofs
Na = 5; % number of CB modes
xi = 0.5e-6; 



%% Newton Method parameters

epsx = 1e-6;
epsf = 1e-6;
maxiter = 100;

%% Coulomb friction

kn = 1e9 ./ (Nx/4); % .* 1e5;
% kn = 3.3e7 ./ (Nx/4); % .* 1e5;
xn0 = 0;
mu = [0.1; 0.1]; % .* 1e5;
kt = [1e6; 1e6] ./ (Nx/4); % .* 1e5;
nloop = 2;

%% preload initial condition for Newton method

xp0 = 0; % if no value defined here, the default value inside is 0


%% continuation parameters


ds = 5;
maxstep = 20000;
% omega_0 = 0.81; % 
% omega_end = 0.87;






% x0 = 0;

% tx0 = 0; % always 0 when calculating the first point so default inside the code
% tomega0 = 1; % is defined by lambda0 and lambda_end

