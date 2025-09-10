clear
load("fc_row.mat");
N = 2^8;
Nx = 4;
xtest = ones(N, 3 * Nx);

tic;
for i = 1:100000
    [F, w] = gf_mex(xtest(1,:), fc);
end
toc;