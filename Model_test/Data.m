
%% External forces and preoload

H_F_ext = [0,0,1]; % fourier coefficient of f(t)


%% HBM parameters

H = 5;
N = 2^8;
Nx = 16;
Na = 3;
xi = 0.05;



%% Newton Method parameters

epsx = 1e-3;
epsf = 1e-3;
maxiter = 100;

%% Coulomb friction

kn = 65.13496132003571;
xn0 = 0;
mu = 0.1 * [1; 1];
kt = 0.06513496132003573 * [1;1];

% w = [-0.803227302337595;-1.06547598368370]; % the middle point of omega 1.2


%% preload initial contition for Newton method

xp0 = 0; % if no value defined here, the default value inside is 0


%% continuation parameters

ds = 0.5;
maxstep = 10000;


omega_0 = 0.6; 
omega_end = 1.1;

x0 = 0;

