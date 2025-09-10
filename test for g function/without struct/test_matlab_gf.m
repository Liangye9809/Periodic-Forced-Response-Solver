clear
load("fc_column.mat");
N = 2^8;
Nx = 4;
xtest = rand(N, 3 * Nx)';
kn = fc.kn;
xn0 = fc.xn0;
mu = fc.mu;
kt = fc.kt;
w_in = fc.w;

tic;
for i = 1:100000
    [F, w] = gf(xtest(:,1), kn, xn0, mu, kt, w_in);
end
toc;

