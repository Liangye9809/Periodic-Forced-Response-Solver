clear
%Build mass matrix (rotor without cyclic matrices)
Mij=load('Mesh_singleBlade2.mas');
index=Mij(:,1)~=Mij(:,2);
ii=[Mij(:,1);Mij(index,2)]; %symmetric part
jj=[Mij(:,2);Mij(index,1)];
ss=[Mij(:,3);Mij(index,3)];
clear Mij;
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
%%
spy(spdiags(M))
Mp = M+M' - sparse(1:length(M),1:length(M),spdiags(M,0));

%%
a = [1,2,3,4,5]';
b = [2,4]';
moveid = ismember(a,b);
keepid = ~moveid;
neworder = [find(keepid); find(moveid)];