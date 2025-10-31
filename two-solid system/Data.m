%% FEM matrices for test
% get the FEM matrices

% for fine mesh

% tep = pwd;
% cd FEM/
% load('CP4.mat');
% cd(tep)
% 
% FEM.Fe = zeros(size(FEM.Mee, 1), 1);
% Fe_idx = [1661, 1667, 31355, 31367] - 1;
% FEM.Fe(Fe_idx) = 1;
% FEM.Fc = zeros(size(FEM.Mcc, 1), 1);

% for old mesh

% tep = pwd;
% cd FEM/two_blades_old/mesh_files/
% load('CP4.mat');
% cd(tep)
% 
% FEM.Fe = zeros(size(FEM.Mee, 1), 1);
% Fe_idx = [112, 118, 5764, 5800];
% FEM.Fe(Fe_idx) = 1;
% FEM.Fc = zeros(size(FEM.Mcc, 1), 1);

% from old input
% tep = pwd;
% cd FEM/
% load('DATA.mat');
% FEM = DATA;
% clear DATA
% cd(tep)
% 
% FEM.Pe = - FEM.Pe;
%% External  forces

H_F_ext = [0,0,1]; % fourier coefficient of f(t)


%% HBM parameters

H = 7; % number of harmonics assumption
N = 2^8; % number of time points per force cycle
Nx = 4; % number of contact points, means having 4 * 3 = 12 dofs
Na = 5; % number of CB modes
xi = 1e-6; 

% calculate the preload forces by predisplacement in normal direction 1 mm

% Xcn0 = [];
% X00 = [0,0,0.001]';
% for i = 1:Nx
%     Xcn0 = [Xcn0; X00];
% end
% FEM.Pe = FEM.Kec * Xcn0;
% FEM.Pc = FEM.Kcc * Xcn0;

%% Newton Method parameters

epsx = 1e-3;
epsf = 1e-3;
maxiter = 100;

%% Coulomb friction

kn = 1e9; % .* 1e5;
xn0 = 0;
mu = [0.1;0.1]; % .* 1e5;
kt = [1e6;1e6]; % .* 1e5;


%% preload initial condition for Newton method

xp0 = 0; % if no value defined here, the default value inside is 0


%% continuation parameters


ds = 0.5;
maxstep = 20000;
omega_0 = 0.81; % 
omega_end = 0.87;
% omega_0 = 4115;
% omega_end = 4318;




% x0 = 0;

% tx0 = 0; % always 0 when calculating the first point so default inside the code
% tomega0 = 1; % is defined by lambda0 and lambda_end

