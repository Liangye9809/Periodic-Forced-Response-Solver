function MSPARSE = read_sparse_ansys(filename)
% READ_SPARSE_ANSYS Reads ANSYS MSPARSE text data into MATLAB sparse matrix.
% It works even when multiple "[i, j]: val" entries appear on one line.

    % Read entire file into a single string
    txt = fileread(filename);

    % Regex pattern capturing:
    %   group 1 → row index
    %   group 2 → col index
    %   group 3 → value
    pattern = '\[\s*(\d+)\s*,\s*(\d+)\s*\]\s*:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)';

    % Extract ALL matches (across lines or space-separated)
    tokens = regexp(txt, pattern, 'tokens');

    if isempty(tokens)
        error("No sparse entries found. Check file formatting.");
    end

    N = numel(tokens);

    % Preallocate triplet arrays
    i = zeros(N,1);
    j = zeros(N,1);
    v = zeros(N,1);

    % Convert text tokens into numbers
    for k = 1:N
        i(k) = str2double(tokens{k}{1});
        j(k) = str2double(tokens{k}{2});
        v(k) = str2double(tokens{k}{3});
    end

    % Determine matrix dimensions
    nRows = max(i);
    nCols = max(j);

    % Build sparse matrix
    MSPARSE = sparse(i, j, v, nRows, nCols);

    fprintf("Loaded %d×%d sparse matrix with %d nonzero entries.\n", ...
        nRows, nCols, nnz(MSPARSE));
end
