clear
%Build mass matrix (rotor without cyclic matrices)
Mij1=load('sinlgeTex/Mesh_Tetra_Single_Linear.mas');
index=Mij1(:,1)~=Mij1(:,2);
ii=[Mij1(:,1);Mij1(index,2)]; %symmetric part
jj=[Mij1(:,2);Mij1(index,1)];
ss=[Mij1(:,3);Mij1(index,3)];
clear Mij1;
M=sparse(ii,jj,ss);
%Build stiffness matrix (rotor without cyclic matrices)
Kij=load('Mesh_singleBlade2.sti');
index=Kij(:,1)~=Kij(:,2);
ii=[Kij(:,1);Kij(index,2)];
jj=[Kij(:,2);Kij(index,1)];
ss=[Kij(:,3);Kij(index,3)];
clear Kij;
K=sparse(ii,jj,ss);
clear ii jj ss;
[Phi, D] = eigs(K, M, 40, 'smallestabs');
e = diag(D);
%% Linear
clear
Mij = load('singleTex/Mesh_Tetra_Single_Linear.mas');
M = sparse(Mij(:, 1), Mij(:, 2), Mij(:, 3));
M = M + M' - sparse(1:length(M), 1:length(M), spdiags(M, 0));
eM = eig(M);
Mfull = full(M);

Kij = load('singleTex/Mesh_Tetra_Single_Linear.sti');
K = sparse(Kij(:, 1), Kij(:, 2), Kij(:, 3));
K = K + K' - sparse(1:length(K), 1:length(K), spdiags(K, 0));
eK = eig(K);
Kfull = full(K);

e = eig(Kfull, Mfull);
%% Quadratic
clear
Mij = load('singleTex/Mesh_4Element_An.mas');
M = sparse(Mij(:, 1), Mij(:, 2), Mij(:, 3));
M = M + M' - sparse(1:length(M), 1:length(M), spdiags(M, 0));
eM = eig(M);
Mfull = full(M);

Kij = load('singleTex/Mesh_4Element_An.sti');
K = sparse(Kij(:, 1), Kij(:, 2), Kij(:, 3));
K = K + K' - sparse(1:length(K), 1:length(K), spdiags(K, 0));
eK = eig(K);
Kfull = full(K);

e = eig(Kfull, Mfull);

%% Linear C3D8
clear
Mij = load('test_box/Mesh_single_element.mas');
M = sparse(Mij(:, 1), Mij(:, 2), Mij(:, 3));
M = M + M' - sparse(1:length(M), 1:length(M), spdiags(M, 0));
eM = eig(M);
Mfull = full(M);

Kij = load('test_box/Mesh_single_element.sti');
K = sparse(Kij(:, 1), Kij(:, 2), Kij(:, 3));
K = K + K' - sparse(1:length(K), 1:length(K), spdiags(K, 0));
eK = eig(K);
Kfull = full(K);

e = eig(Kfull, Mfull);
%% Quadratic C3D20
clear
Mij = load('test_box/Mesh_single_element_Qua.mas');
M = sparse(Mij(:, 1), Mij(:, 2), Mij(:, 3));
M = M + M' - sparse(1:length(M), 1:length(M), spdiags(M, 0));
eM = eig(M);
Mfull = full(M);

Kij = load('test_box/Mesh_single_element_Qua.sti');
K = sparse(Kij(:, 1), Kij(:, 2), Kij(:, 3));
K = K + K' - sparse(1:length(K), 1:length(K), spdiags(K, 0));
eK = eig(K);
Kfull = full(K);

e = eig(Kfull, Mfull);
%%
a = [1,2,3,4,5]';
b = [2,4]';
moveid = ismember(a,b);
keepid = ~moveid;
neworder = [find(keepid); find(moveid)];