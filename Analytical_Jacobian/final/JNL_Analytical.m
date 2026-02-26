function [JNL, JNLt, Ft, wt, Mft] = JNL_Analytical(X, p)
    
    E = p.HBM.E;
    EH = p.HBM.EH;
    Nx = p.HBM.Nx;
    N = p.HBM.N;
    H = p.HBM.H;

    kn = p.fc.kn;
    kt = p.fc.kt;
    xn0 = p.fc.xn0;
    mu = p.fc.mu;
    w_in = p.fc.w;
    nloop = p.fc.nloop;

    xt = zeors(N, 3 * Nx);
    dx = zeors(N, 3 * Nx);
    Xi = zeros(2 * H + 1, 1);
    dXi = zeros(2 * H + 1, 1);

    for i = 1:3 * Nx
        Xi = X((2 * H + 1) * (i - 1) + 1:(2 * H + 1) * i);
        xt(:, i) = E * Xi;

        dXi = dXinFourier(Xi, H);
        dx(:, i) = E * dXi;
    end

    [Ft, wt, Mft] = g(xt, kn, xn0, mu, kt, w_in, nloop);
    
  
    JNLt = zeros(3 * Nx * N, (2 * H + 1) * 3 * Nx);
    JNL = zeros((2 * H + 1) * 3 * Nx, (2 * H + 1) * 3 * Nx);

    for i = 1:Nx
        indxt1 = N * 3 * (i - 1) + 1;
        indxt2 = N * 3 * i;
        indxF1 = (2 * H + 1) * 3 * (i - 1) + 1;
        indxF2 = (2 * H + 1) * 3 * i;

        Mfti = Mft(:, i, :);
        dxi = dx(:, 3 * (i - 1) + 1:3 * i);

        [JNLi, JNLti] = JNL_Analytical_perNx(dxi, kn(i), mu(:, i), kt(:, i), H, N, Mfti, EH);

        JNL(indxF1:indxF2, indxF1:indxF2) = JNLi;
        JNLt(indxt1:indxt2, indxF1:indxF2) = JNLti;

    end


end

function dX = dXinFourier(X, H)
    dX = zeros(size(X));
    for i = 1:H
        dX(2 * i, :) =  i .* X(2 * i + 1, :);
        dX(2 * i + 1, :) =  -i .* X(2 * i, :);
    end

end