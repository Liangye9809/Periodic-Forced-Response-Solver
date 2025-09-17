
%% External  forces and preoload

H_F_ext = [0,0,1]; % fourier coefficient of f(t)


%% HBM parameters

H = 1;
N = 2^8;
Nx = 16;
Na = 3;
xi = 0.05;



%% Newton Method parameters

epsx = 1e-3;
epsf = 1e-3;
maxiter = 100;

%% Coulomb friction

kn = 1e9;
xn0 = 0;
mu = [0.1; 0.1];
kt = [1e7; 1e7];


%% preload initial contition for Newton method

xp0 = 0; % if no value defined here, the default value inside is 0


%% continuation parameters


ds = 0.01;
maxstep = 20000;
omega_0 = 0.6;
omega_end = 1.1;




x0 = 0;


