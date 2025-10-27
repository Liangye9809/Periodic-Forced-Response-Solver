%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INPUT %%%%%%%%%
% M_file_1: mass matrix file (.mas) of solid one (for example blade)
% M_file_2: mass matrix file (.mas) of solid two (for example disk)
% K_file_1: stiffness matrix file (.sti) of silid one
% K_file_2: stiffness matrix file (.sti) of silid two
% dof_file_1: dofs (.dof) file of solid one
% dof_file_1: dofs (.dof) file of solid two
% mesh_file_1: mesh (.inf) file of solid one
% mesh_file_2: mesh (.inf) file of solid two
% node_header_1: nodes' set header of solid one wanted to operate
% node_header_2: nodes' set header of solid two wanted to operate
% CP_1: contact nodes' indices of solid one
% CP_2: contact nodes' indices of solid two
%%%%%%%% OUTPUT %%%%%%%%
% M: mass matrix of whole solid
% K: stiffness matrix of whole solid
%%%%%%%
% CP strores the information of contact nodes, first column is the id of
% contact nodes which needed to reorder, second and third columns are the
% nodes to consist the node surface
function [Mee, Mec, Mcc, Kee, Kec, Kcc, dof] = ReorderMK(M_file_1, M_file_2, K_file_1, K_file_2, ...
                            dof_file_1, dof_file_2, mesh_file_1, mesh_file_2, ...
                            node_header_1, node_header_2, CP_1, CP_2)
%% check the number of contact nodes
if(size(CP_1) ~= size(CP_2))
    error('the number of contact points in two surfaces are not match');
end
Nc = size(CP_1, 1); % number of contacts

%% Build mass and stiffness matrices

%Build mass matrix

Mij = load(M_file_1);
index = Mij(:, 1) ~= Mij(:, 2);
ii = [Mij(:, 1); Mij(index, 2)]; % symmetric part
jj = [Mij(:, 2); Mij(index, 1)];
ss = [Mij(:, 3); Mij(index, 3)];
M1 = sparse(ii, jj, ss);

Mij = load(M_file_2);
index = Mij(:, 1) ~= Mij(:, 2);
ii = [Mij(:, 1); Mij(index, 2)]; % symmetric part
jj = [Mij(:, 2); Mij(index, 1)];
ss = [Mij(:, 3); Mij(index, 3)];
M2 = sparse(ii, jj, ss);
clear Mij;

%Build stiffness matrix

Kij = load(K_file_1);
index = Kij(:, 1) ~= Kij(:, 2);
ii = [Kij(:, 1); Kij(index, 2)];
jj = [Kij(:, 2); Kij(index, 1)];
ss = [Kij(:, 3); Kij(index, 3)];
K1 = sparse(ii,jj,ss);


Kij = load(K_file_2);
index = Kij(:, 1) ~= Kij(:, 2);
ii = [Kij(:, 1); Kij(index, 2)];
jj = [Kij(:, 2); Kij(index, 1)];
ss = [Kij(:, 3); Kij(index, 3)];
K2 = sparse(ii,jj,ss);

clear Kij;
clear ii jj ss;

M = blkdiag(M1, M2); clear M1 M2;
K = blkdiag(K1, K2); clear K1 K2;

% read dof file

fid = fopen(dof_file_1,'r');
dof1 = textscan(fid,'%f');
fclose(fid);
dof1 = dof1{1};  % [1000001.1; 1000001.2; ...]

fid = fopen(dof_file_2,'r');
dof2 = textscan(fid,'%f');
fclose(fid);
dof2 = dof2{1};  % [2000001.1; 2000001.2; ...]

dof = [dof1; dof2]; clear dof1 dof2;


%% permute the matrices and dofs

% separate nodes' id and dofs

node_id = floor(dof);
% dof_dir = round((dof - node_id)*10); % 1, 2, 3

CP = [CP_1; CP_2];
idx_front = find(~ismember(node_id, CP(:, 1)));
idx_back = [];
for i = 1:size(CP,1)
    idx_back = [idx_back; find(node_id == CP(i, 1))];
end

% check the contact dof
if length(idx_back) ~= 3 * 2 * Nc
    error('the contacts dofs in two surfaces are fewer than (3 * 2 * Nc)');
end

new_order = [idx_front; idx_back];


% permute the dof and M and K matrices
dof = dof(new_order);
M = M(new_order, new_order);
K = K(new_order, new_order);

% check the M and K are symmetric still
if issymmetric(M) && issymmetric(K)
    disp('contact nodes put the end of dof in [CP_1; CP_2] order');
else
    error('after permutation - put the contact dofs to the end, M or K loss the symmetricity');
end
% contact nodes at the end of dof in CP order

%% change of basis to tangential and normal


% Read the mesh file1
fid = fopen(mesh_file_1, 'r');
mesh1 = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
mesh1 = mesh1{1};

% extract all the nodes coordinate in *NODE section
node_table_1 = [];
in_node = false;
for i = 1:length(mesh1)
    line = strtrim(mesh1{i});
    
    % Detect section headers
    if startsWith(line, node_header_1, 'IgnoreCase', true) % ignore letter case (uppercase/lowercase)
        in_node = true;
        continue;

    elseif startsWith(line, '*') % Other headers
        in_node = false;
        continue;
    end
    
    % Process node lines
    if in_node && ~isempty(line)
        nums = sscanf(line, '%f,');
        node_table_1 = [node_table_1; nums'];
    end

end
clear mesh1

% Read the mesh file2
fid = fopen(mesh_file_2, 'r');
mesh2 = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
mesh2 = mesh2{1};

node_table_2 = [];
in_node = false;
for i = 1:length(mesh2)
    line = strtrim(mesh2{i});
    
    % Detect section headers
    if startsWith(line, node_header_2, 'IgnoreCase', true) % ignore letter case (uppercase/lowercase)
        in_node = true;
        continue;

    elseif startsWith(line, '*') % Other headers
        in_node = false;
        continue;
    end
    
    % Process node lines
    if in_node && ~isempty(line)
        nums = sscanf(line, '%f,');
        node_table_2 = [node_table_2; nums'];
    end

end
clear mesh2

% find the contact coordinates


coords_1 = zeros(6, 3, Nc);
coords_2 = zeros(6, 3, Nc);

for i = 1:Nc
    for j = 1:3
        % surface 1
        id = CP_1(i, j);
        idx = find(node_table_1(:,1) == id);
        if ~isempty(idx)
            coords_1(j, :, i) = node_table_1(idx, 2:4);
        else
            warning('Node %d not found in file.', id);
        end

        % surface 2
        id = CP_2(i, j);
        idx = find(node_table_2(:,1) == id);
        if ~isempty(idx)
            coords_2(j, :, i) = node_table_2(idx, 2:4);
        else
            warning('Node %d not found in file.', id);
        end
    end
end

% build the nodes' pairs
for i = 1:Nc

    % check if the nodes' pairs match
    for j = 1:3
        d = norm(coords_1(j, : ,i) - coords_2(j, : ,i));
        if d > 1e-8
            error('node pair %d does not match, the distance is %f', i, d);
        end
    end
    % t1 direction
    T1_1 = coords_1(2, :, i) - coords_1(1, :, i);
    t1_1 = T1_1 / norm(T1_1);

    T1_2 = coords_2(2, :, i) - coords_2(1, :, i);
    t1_2 = T1_2 / norm(T1_2);

    % t2 direction
    T2_1 = coords_1(3, :, i) - coords_1(1, :, i);
    t2_1 = T2_1 / norm(T2_1);

    T2_2 = coords_2(3, :, i) - coords_2(1, :, i);
    t2_2 = T2_2 / norm(T2_2);

    % normal direction
    N_1 = cross(t1_1, t2_1);
    n_1 = N_1 / norm(N_1);
    
    N_2 = cross(t1_2, t2_2);
    n_2 = N_2 / norm(N_2);

    coords_1(4:6, :, i) = [t1_1; t2_1; n_1];
    coords_2(4:6, :, i) = [t1_2; t2_2; n_2];

end



%%%%%%%%%%%%%%% full matrix %%%%%%%%%%%%%%%%%
% N_M = size(M, 1);
% Q = eye(N_M);
% for i = 1:Nc
%     idx1 = N_M + 1 - (2 * Nc - i + 1) * 3;
%     idx2 = N_M + 1 - (Nc - i + 1) * 3;
%     T = coords_1(4:6, :, i);
%     Q(idx1:idx1 + 2, idx1:idx1 + 2) = T';
%     Q(idx2:idx2 + 2, idx2:idx2 + 2) = T';
% end
% 
% M = Q' * M * Q;
% K = Q' * K * Q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% block matrix %%%%%%%%%%%%%%%%%
N_Q = Nc * 2 * 3;
Q = zeros(N_Q);

M12 = M(1:end - N_Q, end - N_Q + 1:end);
M22 = M(end - N_Q + 1:end, end - N_Q + 1:end);

K12 = K(1:end - N_Q, end - N_Q + 1:end);
K22 = K(end - N_Q + 1:end, end - N_Q + 1:end);

% build transformation matrix Q
for i = 1:Nc
    idx1 = (i - 1) * 3 + 1;
    idx2 = (i - 1) * 3 + 1 + Nc * 3;
    T = coords_1(4:6, :, i);
    Q(idx1:idx1 + 2, idx1:idx1 + 2) = T';
    Q(idx2:idx2 + 2, idx2:idx2 + 2) = T';
end

% build permutation matrix P - build contact nodes pairs
I = eye(Nc * 3);
P = [I, I; I, -I];

M12 = M12 * Q * P;
M22 = P' * Q' * M22 * Q * P;

K12 = K12 * Q * P;
K22 = P' * Q' * K22 * Q * P;

M(1:end - N_Q, end - N_Q + 1:end) = M12;
M(end - N_Q + 1:end, end - N_Q + 1:end) = M22;
M(end - N_Q + 1:end, 1:end - N_Q) = M12';

K(1:end - N_Q, end - N_Q + 1:end) = K12;
K(end - N_Q + 1:end, end - N_Q + 1:end) = K22;
K(end - N_Q + 1:end, 1:end - N_Q) = K12';

% check the M and K are symmetric still
if issymmetric(M) && issymmetric(K)
    disp('contact basis on tangential and normal in pairs');
else
    error('after transformation - change contact basis to tangential and normal, M or K loss the symmetricity');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% split M and K
Mee = M(1:end - Nc * 3, 1:end - Nc * 3);
Mec = M(1:end - Nc * 3, end - Nc * 3 + 1:end);
Mcc = M(end - Nc * 3 + 1:end, end - Nc * 3 + 1:end);

Kee = K(1:end - Nc * 3, 1:end - Nc * 3);
Kec = K(1:end - Nc * 3, end - Nc * 3 + 1:end);
Kcc = K(end - Nc * 3 + 1:end, end - Nc * 3 + 1:end);

end

