%% for blade 1
B1Nodes = readmatrix('Sub-mesh_ContactFace_1.dat');
A = B1Nodes; % make a copy
B1M = zeros(65, 65); % final form of mesh matrix
for i = 1:64
    ymaxi = max(A(:,3));
    idx = A(:,3) == ymaxi;
    Gi = A(idx, :); % get the i-th row of mesh nodes
    A = A(~idx, :); % delete the i-th row of mesh nodes
    [~, sortedIdx] = sort(Gi(:, 2)); % sort the x-coordinate
    Gi = Gi(sortedIdx, :);
    B1M(i, :) = Gi(:, 1)';
end
[~, sortedIdx] = sort(A(:, 2)); % sort the x-coordinate
A = A(sortedIdx, :);
B1M(65, :) = A(:, 1)';

%% for blade 2
B2Nodes = readmatrix('Sub-mesh_ContactFace_2.dat');
A = B2Nodes; % make a copy
B2M = zeros(65, 65); % final form of mesh matrix
for i = 1:64
    ymaxi = max(A(:,3));
    idx = A(:,3) == ymaxi;
    Gi = A(idx, :); % get the i-th row of mesh nodes
    A = A(~idx, :); % delete the i-th row of mesh nodes
    [~, sortedIdx] = sort(Gi(:, 2)); % sort the x-coordinate
    Gi = Gi(sortedIdx, :);
    B2M(i, :) = Gi(:, 1)';
end
[~, sortedIdx] = sort(A(:, 2)); % sort the x-coordinate
A = A(sortedIdx, :);
B2M(65, :) = A(:, 1)';