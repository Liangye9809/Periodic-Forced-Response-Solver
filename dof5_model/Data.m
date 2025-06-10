%% FEM matrices for test
FEM.Mee = [1,0;0,1];
FEM.Mec = [0,0,0;0,0,0];
FEM.Mcc = [1,0,0;0,1,0;0,0,1];
FEM.Kee = [10,-10;-10,20];
FEM.Kec = [0,0,0;-10,0,0];
FEM.Kcc = [20,-10,0;-10,20,-10;0,-10,20]; 
FEM.Fe = [1,1]'; % external forces amplitude in elastic part, Ne * 1 vector
FEM.Fc = [0,0,0]'; % external forces amplitude in contact part, Nc * 1 vector
FEM.Pe = [1;0.1];
FEM.Pc = [0.1;0.1;0.1];

%% CB matrices for test
% Phi = [-0.850650808352040,-0.525731112119134;-0.525731112119134,0.850650808352040];
% Psi = [1,0,0;1,0,0];
% Kaa = [3.819660112501051;26.180339887498950];
% Kxx = [10,-10,0;-10,20,-10;0,-10,20];
% Max = [-1.376381920471173,0,0;0.324919696232906,0,0];
% Mxx = [3,0,0;0,1,0;0,0,1];
% Fa = [-1.376381920471173;0.324919696232906];
% Fx = [2;0;0];
% Rx = [-1.2;-0.1;-0.1];
% xe0 = [0.21;0.11];

%% External forces and preoload

H_F_ext = [0,0,1]; % fourier coefficient of f(t)


%% HBM parameters

H = 10;
N = 2^5;
Nx = 1; % contact points
Na = 2;
xi = 0.05;

%% Newton Method parameters

epsx = 1e-5;
epsf = 1e-5;
maxiter = 50;

%% Coulomb friction

kn = 1;
xn0 = -3.5;
mu = [1;1];
kt = [1;1];

w = [-0.803227302337595;-1.06547598368370]; % the middle point of omega 1.2
%% preload initial contition for Newton method

xp0 = 1e-5; % if no value defined here, the default value inside is 0


%% continuation parameters


ds = 0.05;
maxstep = 10000;
omega_0 = 1.2;
omega_end = 1.3;

cd ./data/
load('x_omega1.2_backwards.mat');
load('x_omega1.2_forwards.mat');
x0 = (x_backwards + x_forwards)/2;
cd ../
% tx0 = 0; % always 0 when calculating the first point so default inside the code
% tomega0 = 1; % is defined by lambda0 and lambda_end

