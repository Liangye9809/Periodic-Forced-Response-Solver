%% FEM matrices for test
% FEM.Mee = [1,0;0,1];
% FEM.Mec = [0,0,0;0,0,0];
% FEM.Mcc = [1,0,0;0,1,0;0,0,1];
% FEM.Kee = [10,-10;-10,20];
% FEM.Kec = [0,0,0;-10,0,0];
% FEM.Kcc = [20,-10,0;-10,20,-10;0,-10,20]; 
% FEM.Fe = [1,1]'; % external forces amplitude in elastic part, Ne * 1 vector
% FEM.Fc = [0,0,0]'; % external forces amplitude in contact part, Nc * 1 vector
% FEM.Pe = [1;0.1];
% FEM.Pc = [0.1;0.1;0.1];

%% CB matrices for test
% CB.CBmods.Phi = [-0.850650808352040,-0.525731112119134;-0.525731112119134,0.850650808352040];
% CB.CBmods.Psi = [1,0,0;1,0,0];
% CB.CBmods.Kaa = [3.819660112501051;26.180339887498950];
% CB.CB_MK.Kxx = [10,-10,0;-10,20,-10;0,-10,20];
% CB.CB_MK.Max = [-1.376381920471173,0,0;0.324919696232906,0,0];
% CB.CB_MK.Mxx = [3,0,0;0,1,0;0,0,1];
% CB.CB_F.Fa = [-1.376381920471173;0.324919696232906];
% CB.CB_F.Fx = [2;0;0];
% CB.Rx = [-1.2;-0.1;-0.1];
% CB.xe0 = [0.21;0.11];

%% External  forces

H_F_ext = [0,0,1]; % fourier coefficient of f(t)


%% HBM parameters

H = 10; % number of harmonics assumption
N = 2^8; % number of time points per force cycle
Nx = 4; % number of contact points, means having 4 * 3 = 12 dofs
Na = 5; % number of CB modes
xi = 1e-6; 



%% Newton Method parameters

epsx = 1e-3;
epsf = 1e-3;
maxiter = 100;

%% Coulomb friction

kn = 1e9;
xn0 = 0;
mu = [0.1;0.1];
kt = [1e6;1e6];


%% preload initial condition for Newton method

xp0 = 0; % if no value defined here, the default value inside is 0


%% continuation parameters


ds = 0.5;
maxstep = 20000;
omega_0 = 0.81; % 
omega_end = 0.85;
% omega_0 = 4115;
% omega_end = 4318;




% x0 = 0;

% tx0 = 0; % always 0 when calculating the first point so default inside the code
% tomega0 = 1; % is defined by lambda0 and lambda_end

