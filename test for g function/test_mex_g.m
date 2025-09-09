clear
load("ptest.mat");
N = 2^8;
Nx = 4;
xtest = zeros(N, 3 * Nx);
fc.kn = ptest.fc.kn;
fc.xn0 = ptest.fc.xn0;
fc.mu = ptest.fc.mu;
fc.kt = ptest.fc.kt;
fc.w = ptest.fc.w;

profile on
tic;
for i = 1:200
    F = g_mex(xtest, fc);
end
toc;
profile off
p = profile('info');
profview(0,p);