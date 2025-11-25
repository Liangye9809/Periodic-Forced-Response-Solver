clear
clc

M = read_sparse_ansys('M.mass');
Mp = M + M' - sparse(1:length(M), 1:length(M), spdiags(M, 0));
M_full = full(Mp);

K = read_sparse_ansys('K.sti');
Kp = K + K' - sparse(1:length(K), 1:length(K), spdiags(K, 0));
K_full = full(Kp);


[Phi, D] = eigs(K_full, M_full, 20, 'smallestabs');
e = diag(D);

[Phi_M, D_M] = eigs(M_full, 20, 'smallestabs');
e_M = diag(D_M);

[Phi_K, D_K] = eig(K_full);
e_K = diag(D_K);

%%
clear
clc
[M, Mv] = hb_to_msm('M_HW.mas');
Mp = M + M' - sparse(1:length(M), 1:length(M), spdiags(M, 0));
M_full = full(Mp);

[K, Kv] = hb_to_msm('K_HW.sti');
Kp = K + K' - sparse(1:length(K), 1:length(K), spdiags(K, 0));
K_full = full(Kp);

[Phi, D] = eigs(K_full, M_full, 20, 'smallestabs');
e = diag(D);
e_Hz = sqrt(e) / (2 * pi);

[Phi_M, D_M] = eigs(M_full, 20, 'smallestabs');
e_M = diag(D_M);

[Phi_K, D_K] = eig(K_full);
e_K = diag(D_K);