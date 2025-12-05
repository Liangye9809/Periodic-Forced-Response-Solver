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
[M, Mv] = hb_to_msm('M.mas');
Mp = M + M' - sparse(1:length(M), 1:length(M), spdiags(M, 0));
M_full = full(Mp);



[K, Kv] = hb_to_msm('K.sti');
Kp = K + K' - sparse(1:length(K), 1:length(K), spdiags(K, 0));
K_full = full(Kp);


% [Phi, D] = eigs(K_full, M_full, 50, 'smallestabs');
% e = diag(D);
% e_Hz = sqrt(e) / (2 * pi);


[Phi_M, D_M] = eig(M_full);
e_M = diag(D_M);


[Phi_K, D_K] = eig(K_full);
e_K = diag(D_K);


[M_keyopt2, Mv_keyopt2] = hb_to_msm('M_keyopt2.mas');
Mp_keyopt2 = M_keyopt2 + M_keyopt2' - sparse(1:length(M_keyopt2), 1:length(M_keyopt2), spdiags(M_keyopt2, 0));
M_full_keyopt2 = full(Mp_keyopt2);


[K_keyopt2, Kv_keyopt2] = hb_to_msm('K_keyopt2.sti');
Kp_keyopt2 = K_keyopt2 + K_keyopt2' - sparse(1:length(K_keyopt2), 1:length(K_keyopt2), spdiags(K_keyopt2, 0));
K_full_keyopt2 = full(Kp_keyopt2);

[Phi_M_keyopt2, D_M_keyopt2] = eig(M_full_keyopt2);
e_M_keyopt2 = diag(D_M_keyopt2);

[Phi_K_keyopt2, D_K_keyopt2] = eig(K_full_keyopt2);
e_K_keyopt2 = diag(D_K_keyopt2);