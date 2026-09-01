function F = Getgx(xp, fc)
    % [F, ~] = g_mex(xp, fc);
       kn = fc.kn;
      xn0 = fc.xn0;
       mu = fc.mu;
       kt = fc.kt;
     w_in = fc.w;
    nloop = fc.nloop;
    [Fi, ~, ~] = g(xp, kn, xn0, mu, kt, w_in, nloop);
    F = Fi(end, :);
end

