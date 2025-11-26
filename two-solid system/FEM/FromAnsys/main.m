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

[M_lump, Mv_lump] = hb_to_msm('M_lumpm.mas');
Mp_lump = M_lump + M_lump' - sparse(1:length(M_lump), 1:length(M_lump), spdiags(M_lump, 0));
M_full_lump = full(Mp_lump);


[K, Kv] = hb_to_msm('K_HW.sti');
Kp = K + K' - sparse(1:length(K), 1:length(K), spdiags(K, 0));
K_full = full(Kp);

[K_lump, Kv_lump] = hb_to_msm('K_lumpm.sti');
Kp_lump = K_lump + K_lump' - sparse(1:length(K_lump), 1:length(K_lump), spdiags(K_lump, 0));
K_full_lump = full(Kp_lump);

[Phi, D] = eigs(K_full, M_full, 50, 'smallestabs');
e = diag(D);
e_Hz = sqrt(e) / (2 * pi);

[Phi_lump, D_lump] = eigs(K_full, M_full_lump, 50, 'smallestabs');
e_lump = diag(D_lump);
e_lump_Hz = sqrt(e_lump) / (2 * pi);

[Phi_M, D_M] = eig(M_full);
e_M = diag(D_M);

[Phi_M_lump, D_M_lump] = eig(M_full_lump);
e_M_lump = diag(D_M_lump);

[Phi_K, D_K] = eig(K_full);
e_K = diag(D_K);

[Phi_K_lump, D_K_lump] = eig(K_full_lump);
e_K_lump = diag(D_K_lump);