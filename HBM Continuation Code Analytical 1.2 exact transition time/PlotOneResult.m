function [ft, wt, flag, xct, dxct, t, XcM] = PlotOneResult(X, w_in, Omega, params)
    Na = params.func.HBM.Na;
    Nx = params.func.HBM.Nx;
    H = params.func.HBM.H;
    N = params.func.HBM.N;
    E = params.func.HBM.E;
    EH = params.func.HBM.EH;
    xp = params.func.static.preload.xp;
    gxp = params.func.static.preload.gxp;
    xn0 = params.func.fc.xn0;
    kt = params.func.fc.kt;
    kn = params.func.fc.kn;
    mu = params.func.fc.mu;
    nloop = params.func.fc.nloop;
    Xc = X((2 * H + 1) * Na + 1:end);
    xct = Fourier_to_Time(Xc, H, Nx, E);
    params.func.fc.w = w_in;
    [FUN, w_poss, ~, flag, S] = HBMFUNC(X, xct + xp', Omega, params.func);
    omega_poss = Omega;
    if norm(FUN) > params.Newton.epsf
        warning('norm(FUN_poss) > params.Newton.epsf');
        params.cont.ds = 0;
        params.cont.omega_0 = Omega;
        params.cont.step = 100001;
        [Xc, omega_poss, ~, ~, w_poss, k_poss, ~, ~] = cont_step(@HBMFUNC, @HBMJACOB, @HBMJOmega, params);
        params.func.fc.w = w_poss;
    end
    xct = Fourier_to_Time(Xc, H, Nx, E);
    dxct = get_dxt(Xc, E, Nx);
    [ft, wt, flag] = g(xct + xp', kn, xn0, mu, kt, w_poss, nloop, dxct);
    % ft = ft - gxp';
    T = 2 * pi / omega_poss;
    dt = T / N;
    t = [0:dt:(T - dt)]';
    XcM = zeros((2 * H + 1), Nx * 3);
    for i = 1:Nx * 3
        XcM(:, i) = Xc((2 * H + 1) * (i - 1) + 1:(2 * H + 1) * i);
    end
end