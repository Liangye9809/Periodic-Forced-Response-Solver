input_file = 'BladeMesh_2.inp';
output_file = 'BladeMesh_2_shifted.inp';
first_section_header = '*NODE, NSET=NALL';
offset = 2e6;
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
