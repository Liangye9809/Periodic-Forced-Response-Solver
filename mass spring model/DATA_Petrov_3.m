CB = [];
CB.CB_MK.Max = zeros(1, 3);
CB.CB_MK.Mxa = zeros(3, 1);
CB.CB_MK.Mxx = diag([1, 1, 1]);

CB.CB_MK.Kaa = 1;
CB.CB_MK.Kxx = diag([135, 40, 40]);

CB.CB_F.Fa = 0;
CB.CB_F.Fx = [100, 0, 100]';

CB.CBmods = [];

H_F_ext = [0, 0, 1];


H = 10; 
N = 2^8; 
Nx = 1; 
Na = 1; 
xi = 0.01; 

epsx = 1e-9;
epsf = 1e-9;
maxiter = 100;

kn = 120; 
xn0 = 0;
mu = 15 * [1; 1]; 
kt = 1 * [1; 1]; 
nloop = 2;

ds = 0.02;
maxstep = 63580;

xp0 = 0;
Rx = [0, 0, -800]';
% xp = [0, 0, 7.5]';
% gxp = [0, 0, 300]';
xe0 = 0;