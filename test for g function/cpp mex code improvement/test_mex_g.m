clear
load("fc_row.mat");
N = 2^8;
Nx = 4;
xtest = rand(N, 3 * Nx);

tic;
for i = 1:1000
    [F, w] = g_mex(xtest, fc);
end
toc;
