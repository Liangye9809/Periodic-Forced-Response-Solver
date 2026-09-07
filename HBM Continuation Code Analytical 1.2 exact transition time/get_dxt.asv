function dxt = get_dxt(X, E, Nx)
    
    n = size(E, 2); % 2H+1
    a = 3 * Nx; % number of DOF
    X_ = zeros(n, a);
    for i = 1:a
        r1 = n * (i - 1) + 1;
        r2 = n * i;
        X_(:,i) = X(r1:r2); % reorder in dofs in column
    end
    dX = zeros(size(X_));
    for h = 1:(n - 1) / 2
        dX(2 * h, :) = h * X_(2 * h + 1, :);
        dX(2 * h + 1, :) = -h * X_(2 * h, :);
    end
    dxt = E * dX; % dx/dtheta; Omega cancels in the slip-exit interpolation
end