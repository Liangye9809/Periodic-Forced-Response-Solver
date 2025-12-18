input_file = 'MeshWithPreload_32x32_1.inp';
output_file = 'MeshWithPreload_32x32_1_shifted.inp';
first_section_header = '*NODE, NSET=NALL';
offset = 1e6;
RenumberMesh(input_file, output_file, first_section_header, offset);
%%
line = '10, 0.00215423, 0.1123, 0.0001';
nums = sscanf(line, '%f,');
%% contact nodes 
% in blade 1
Nc1 = [153, 3915, 3952;
       148, 3900, 3937;
       105, 3763, 3800;
       100, 3748, 3785] + 1e6;
% in balde 2
Nc2 = [499, 3427, 3467;
       504, 3442, 3482;
       451, 3275, 3315;
       456, 3290, 3330] + 2e6;
%%
M1 = rand(3);
M2 = 10 * rand(4);
blkdiag(M1, M2)
%%
move_idx = [0,1,1,1,1,0,0]';
keep_idx = ~move_idx;
find(move_idx);

%% new 2 solid mesh
clear
M_file_1 = 'Blade1Anlysis.mas';
M_file_2 = 'Blade2Anlysis.mas';
K_file_1 = 'Blade1Anlysis.sti';
K_file_2 = 'Blade2Anlysis.sti';
dof_file_1 = 'Blade1Anlysis.dof';
dof_file_2 = 'Blade2Anlysis.dof';
mesh_file_1 = 'BladeMesh_1_shifted.inp';
mesh_file_2 = 'BladeMesh_2_shifted.inp';
node_header_1 = '*NODE';
node_header_2 = '*NODE';
Nc1 = [153, 3915, 3952;
       148, 3900, 3937;
       105, 3763, 3800;
       100, 3748, 3785; %] + 1e6; % CP4
       179, 3997, 4023;
       170, 3970, 4005;
        83, 3693, 3730;
        74, 3666, 3703; %] + 1e6; % CP8
       131, 3845, 3882;
       122, 3818, 3855;
       175, 3985, 4015;
        78, 3678, 3715; %] + 1e6; % CP12
        61, 3924, 1772;
        70, 3888, 1782;
        57, 3772, 1768;
        66, 3736, 1778; %] + 1e6; % CP16
        39, 1708, 4019;
        34, 1703, 4009;
        51, 1761, 3686;
        46, 1756, 3671; %] + 1e6; % CP20
       155, 3921, 3958;
       146, 3894, 3931;
       107, 3769, 3806;
        98, 3742, 3779; %] + 1e6; % CP24
         6, 1711, 1774;
         4, 1699, 1784;
         2, 1764, 1765;
         1, 1752, 1775; %] + 1e6; % CP28
       129, 3839, 3876;
       124, 3824, 3861;
       150, 3906, 3943;
       103, 3757, 3794; %] + 1e6; % CP32
        59, 3848, 1770;
        68, 3812, 1780;
       139, 3871, 3908;
       114, 3792, 3829] + 1e6; % CP36
% in balde 2
Nc2 = [499, 3427, 3467;
       504, 3442, 3482;
       451, 3275, 3315;
       456, 3290, 3330; %] + 2e6; % CP4
       521, 3497, 3535;
       530, 3524, 3553;
       425, 3193, 3233;
       434, 3220, 3260; %] + 2e6; % CP8
       473, 3345, 3385;
       482, 3372, 3412;
       525, 3509, 3543;
       430, 3208, 3248; %] + 2e6; % CP12
       421, 3418, 1597;
       388, 3454, 1365;
       417, 3266, 1593;
       384, 3302, 1361; %] + 2e6; % CP16
       406, 1487, 3539;
       411, 1492, 3549;
       394, 1375, 3201;
       399, 1380, 3216; %] + 2e6; % CP20
       497, 3421, 3461;
       506, 3448, 3488;
       449, 3269, 3309;
       458, 3296, 3336; %] + 2e6; % CP24
        13, 1484, 1599;
        30, 1496, 1367;
        24, 1372, 1590;
        27, 1384, 1358; %] + 2e6; % CP28
       475, 3351, 3391;
       480, 3366, 3406;
       502, 3436, 3476;
       453, 3281, 3321; %] + 2e6; % CP32
       419, 3342, 1595;
       386, 3378, 1363;
       489, 3395, 3435;
       466, 3322, 3362] + 2e6;
% [FEM.Mee, FEM.Mec, FEM.Mcc, FEM.Kee, FEM.Kec, FEM.Kcc, dof] = ReorderMK(M_file_1, M_file_2, ...
%                    K_file_1, K_file_2, ...
%                    dof_file_1, dof_file_2, ...
%                    mesh_file_1, mesh_file_2, ...
%                    node_header_1, node_header_2, ...
%                    Nc1, Nc2);
% [Phi, D] = eigs(FEM.Kee, FEM.Mee, 20, 'smallestabs');
% e = diag(D);
% save CP36 FEM
% save full_stuck_modes_CP36 e
% Fe_dof = [1000572.2, 1000574.2, 2000596.2, 2000592.2]';
% idx_Fe = find(ismember(dof, Fe_dof));

for i = 4:4:36
    Nx_1 = Nc1(1:i, :);
    Nx_2 = Nc2(1:i, :);
    [FEM.Mee, FEM.Mec, FEM.Mcc, FEM.Kee, FEM.Kec, FEM.Kcc, dof] = ReorderMK(M_file_1, M_file_2, ...
                       K_file_1, K_file_2, ...
                       dof_file_1, dof_file_2, ...
                       mesh_file_1, mesh_file_2, ...
                       node_header_1, node_header_2, ...
                       Nx_1, Nx_2);
    Pe_dof1 = [263; 262; 14; 809; 199; 208; 209; 15] + 0.3 + 1e6;
    Pe_dof2 = [121; 120; 29; 236; 119; 117; 118; 28] + 0.3 + 2e6;
    FEM.idx_Pe1 = find(ismember(dof, Pe_dof1));
    FEM.idx_Pe2 = find(ismember(dof, Pe_dof2));

    Fe_dof = [1000572.2, 1000574.2, 2000596.2, 2000592.2]';
    FEM.idx_Fe = find(ismember(dof, Fe_dof));
    % save CP36 FEM
    save('Mesh10x13_CP' + string(i) + '.mat', 'FEM');
end
%%
M_file_1 = 'StructralMesh_8X8_1_shifted.mas';
M_file_2 = 'StructralMesh_8X8_2_shifted.mas';
K_file_1 = 'StructralMesh_8X8_1_shifted.sti';
K_file_2 = 'StructralMesh_8X8_2_shifted.sti';
dof_file_1 = 'StructralMesh_8X8_1_shifted.dof';
dof_file_2 = 'StructralMesh_8X8_2_shifted.dof';
mesh_file_1 = 'StructralMesh_8X8_1_shifted.inp';
mesh_file_2 = 'StructralMesh_8X8_2_shifted.inp';
node_header_1 = '*NODE';
node_header_2 = '*NODE';
Nc1 = [57, 9248, 9249] + 1e6;
       % 20, 7, 26;
       % 10, 26, 1;
       % 26, 14, 18] + 1e6;
% in balde 2
Nc2 = [57, 7881, 7906] + 2e6;
       % 19, 8, 25;
       % 12, 25, 2;
       % 25, 16, 17] + 2e6;
[FEM.Mee, FEM.Mec, FEM.Mcc, FEM.Kee, FEM.Kec, FEM.Kcc, dof] = ReorderMK(M_file_1, M_file_2, ...
                   K_file_1, K_file_2, ...
                   dof_file_1, dof_file_2, ...
                   mesh_file_1, mesh_file_2, ...
                   node_header_1, node_header_2, ...
                   Nc1, Nc2);
[Phi, D] = eigs(FEM.Kee, FEM.Mee, 40, 'smallestabs');
e = diag(D);
save('StructuralCP' + string(1) + '.mat', 'FEM');
%%
I = eye(6);
P = [I, I; I, -I];
P = inv(P);
%%
filename = '/mnt/data/liangye/Downloads/two_blades/inputs/input_CBtest.txt';  % change to your filename
% Read the mesh file1
fid = fopen(filename, 'r');
data = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
data = data{1};

% extract all the nodes coordinate in *NODE section
DATA = struct();
entries = [];
currentName = '';
% Name = false;
for i = 1:length(data)
    line = strtrim(data{i});
    
    % Detect section headers
    if startsWith(line, 'Name:', 'IgnoreCase', true) % ignore letter case (uppercase/lowercase)
        if ~isempty(currentName)
            DATA.(currentName) = entries;
        end
        currentName = strtrim(extractAfter(line, 'Name:'));
        % Name = true;
        entries = [];
        continue;
    end
    
    % Process node lines
    if ~isempty(line)
        nums = sscanf(line, '%f,');
        entries = [entries; nums'];
    end

end
save DATA.mat DATA
% load("DATA.mat");
DATA.Kec = sparse(DATA.Kec(:,1), DATA.Kec(:,2), DATA.Kec(:,3));
DATA.Kcc = sparse(DATA.Kcc(:,1), DATA.Kcc(:,2), DATA.Kcc(:,3));
DATA.Kee = sparse(DATA.Kee(:,1), DATA.Kee(:,2), DATA.Kee(:,3));

DATA.Mec = sparse(DATA.Mec(:,1), DATA.Mec(:,2), DATA.Mec(:,3));
DATA.Mcc = sparse(DATA.Mcc(:,1), DATA.Mcc(:,2), DATA.Mcc(:,3));
DATA.Mee = sparse(DATA.Mee(:,1), DATA.Mee(:,2), DATA.Mee(:,3));

DATA.Fe = sparse(DATA.Fe(:,1), DATA.Fe(:,2), DATA.Fe(:,3));

DATA.Pc = sparse(DATA.Pc(:,1), DATA.Pc(:,2), DATA.Pc(:,3));
DATA.Pe = - Pe;

Pe = - Kec * Xcn0;
Pc = - Kcc * Xcn0;


DATA.Fc = zeros(12,1);
Fe = zeros(size(DATA.Mee,1),1);
Fe(1:5800) = DATA.Fe;
DATA.Fe = sparse(Fe);
DATA.Pe = sparse(DATA.Pe);

Kec = zeros(size(DATA.Mee,1),12);
Mec = zeros(size(DATA.Mee,1),12);
Kec(1:10563,:) = DATA.Kec;
Mec(1:10563,:) = DATA.Mec;
DATA.Kec = sparse(Kec);
DATA.Mec = sparse(Mec);
Xcn0 = [0,0,0.001,0,0,0.001,0,0,0.001,0,0,0.001]';
DATA.Pe = - DATA.Kec * Xcn0;

%% for old mesh
M_file_1 = 'two_blades_old/mesh_files/analysis_left.mas';
M_file_2 = 'two_blades_old/mesh_files/analysis_right.mas';
K_file_1 = 'two_blades_old/mesh_files/analysis_left.sti';
K_file_2 = 'two_blades_old/mesh_files/analysis_right.sti';
dof_file_1 = 'two_blades_old/mesh_files/analysis_left.dof';
dof_file_2 = 'two_blades_old/mesh_files/analysis_right.dof';
mesh_file_1 = 'two_blades_old/mesh_files/blade-left-mesh-clean.inp';
mesh_file_2 = 'two_blades_old/mesh_files/blade-right-mesh-clean.inp';
node_header_1 = '*NODE';
node_header_2 = '*NODE';
Nc1 = [545, 23, 543;
       548, 24, 546;
       529, 19, 527;
       532, 20, 530] + 1e6;
       
% in balde 2
Nc2 = [43, 33, 40;
       47, 36, 44;
       63, 51, 60;
       67, 55, 64] + 2e6;
       
[FEM.Mee, FEM.Mec, FEM.Mcc, FEM.Kee, FEM.Kec, FEM.Kcc, dof] = ReorderMK(M_file_1, M_file_2, ...
                   K_file_1, K_file_2, ...
                   dof_file_1, dof_file_2, ...
                   mesh_file_1, mesh_file_2, ...
                   node_header_1, node_header_2, ...
                   Nc1, Nc2);
[Phi, D] = eigs(FEM.Kee, FEM.Mee, 40, 'smallestabs');
e = diag(D);
save two_blades_old/mesh_files/CP4 FEM

%% mesh 32x32 contact with preload mesh 2x2(9 nodes)
clear
M_file_1 = 'MeshWithPreload_32x32_1_shifted.mas';
M_file_2 = 'MeshWithPreload_32x32_2_shifted.mas';
K_file_1 = 'MeshWithPreload_32x32_1_shifted.sti';
K_file_2 = 'MeshWithPreload_32x32_2_shifted.sti';
dof_file_1 = 'MeshWithPreload_32x32_1_shifted.dof';
dof_file_2 = 'MeshWithPreload_32x32_2_shifted.dof';
mesh_file_1 = 'MeshWithPreload_32x32_1_shifted.inp';
mesh_file_2 = 'MeshWithPreload_32x32_2_shifted.inp';
node_header_1 = '*NODE';
node_header_2 = '*NODE';

% % 1 contact
% Nc1 = [618, 9458, 9552] + 1e6;
% Nc2 = [618, 9465, 9469] + 2e6;

% 4 contact
% Nc1 = [858, 10194, 10288;
%        874, 10242, 10336;
%        362, 8674, 8768;
%        378, 8722, 8816] + 1e6;
% 
% Nc2 = [874, 10249, 10253;
%        858, 10201, 10205;
%        378, 8729, 8733;
%        362, 8681,8685] + 2e6;

% 16 contact
% Nc1 = [978, 10562, 10656;
%        986, 10586, 10680;
%        1002, 10634, 10728;
%        730, 9802, 9896;
%        738, 9826, 9920;
%        746, 9850, 9944;
%        754, 9874, 9968;
%        482, 9042, 9136;
%        490, 9066, 9160;
%        498, 9090, 9184;
%        506, 9114, 9208;
%        234, 8282, 8376;
%        424, 8306, 8400;
%        250, 8330, 8424;
%        258, 8354, 8448] + 1e6;
% 
% Nc2 = [] + 2e6;

% calculate contact pairs
load B1M.mat
load B2M.mat
Nc1 = [];
Nc2 = [];
NcTemp = zeros(1, 3);
for i = 1:5
    n = 2^i;
    m = 2^(i-1);
    
    for j = (64 / n + 1):(64 / m):64
        for k = (64 / n + 1):(64 / m):64
            NcTemp(1) = B1M(j, k);
            NcTemp(2) = B1M(j, k + 1);
            NcTemp(3) = B1M(j - 1, k);
            Nc1 = [Nc1; NcTemp];

            NcTemp(1) = B2M(j, k);
            NcTemp(2) = B2M(j, k + 1);
            NcTemp(3) = B2M(j - 1, k);
            Nc2 = [Nc2; NcTemp];
        end
    end
end

Nc1 = Nc1 + 1e6;
Nc2 = Nc2 + 2e6;
for i = 1:5
    S1 = (4^(i - 1) - 1) / 3 + 1;
    S1 = max(1, S1);
    S2 = (4^i - 1) / 3;
    [FEM.Mee, FEM.Mec, FEM.Mcc, FEM.Kee, FEM.Kec, FEM.Kcc, dof] = ReorderMK(M_file_1, M_file_2, ...
                       K_file_1, K_file_2, ...
                       dof_file_1, dof_file_2, ...
                       mesh_file_1, mesh_file_2, ...
                       node_header_1, node_header_2, ...
                       Nc1(S1:S2,:), Nc2(S1:S2,:));
    
    Pe_dof1 = [1; 2; 3; 4; 5; 6; 7; 8; 9] + 0.3 + 1e6;
    Pe_dof2 = [1; 2; 3; 4; 5; 6; 7; 8; 9] + 0.3 + 2e6;
    FEM.idx_Pe1 = find(ismember(dof, Pe_dof1));
    FEM.idx_Pe2 = find(ismember(dof, Pe_dof2));
    
    Fe_dof = [1001432.2, 1001462.2, 2001719.2, 2001765.2]';
    FEM.idx_Fe = find(ismember(dof, Fe_dof));
    
    % [Phi, D] = eigs(FEM.Kee, FEM.Mee, 20, 'smallestabs');
    % e = diag(D);
    % save full_stuck_modes_CP36 e
    
    % save Mesh32x32_CP4 FEM
    save('Mesh32x32_CP' + string(4^(i-1)) + '.mat', 'FEM');
end

%%
