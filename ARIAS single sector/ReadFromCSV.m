path = pwd;
% cd('D:\study\PHD\data\Frictions\HBM Keep A and X dof'); % for windows at home
% cd('E:\study loves me\UPM\PHD\data\non-linear problem\Frictions\HBM Keep A and X dof'); % for windows laptop
cd('/mnt/data/liangye/Downloads/ARIAS') % for Liunx
filename = 'data.csv';
fid = fopen(filename, 'r');

matrices = struct();  % Store matrices here
currentName = '';
dataBlock = {};

while ~feof(fid)
    line = strtrim(fgetl(fid));

    if startsWith(line, 'Name:')
        % Save previous matrix if exists
        if ~isempty(currentName)
            matrices.(currentName) = cell2mat(dataBlock);
        end

        % Clean and validate field name
        rawName = strrep(line, 'Name:', '');
        currentName = matlab.lang.makeValidName(strtrim(rawName));
        dataBlock = {};
    elseif ~isempty(line)
        row = str2double(strsplit(line, ','));
        dataBlock{end+1,1} = row;
    end
end

% Save the last matrix
if ~isempty(currentName)
    matrices.(currentName) = cell2mat(dataBlock);
end

fclose(fid);

% Access with:
% matrices.Mxx
% matrices.Kxx

cd(path);

CB.CBmods.Phi = matrices.Phi;
CB.CBmods.Psi = matrices.Psi;
CB.CBmods.Kaa = matrices.omega_a2;

CB.CB_MK.Max = matrices.Max;
CB.CB_MK.Mxx = matrices.Mxx;
CB.CB_MK.Kaa = matrices.omega_a2;
CB.CB_MK.Kxx = matrices.Kxx;

CB.CB_F.Fa = matrices.fa;
CB.CB_F.Fx = matrices.fx;

xe0 = matrices.xe0;
Rx = matrices.Rx;

save CB.mat CB
save xe0.mat xe0
save Rx.mat Rx