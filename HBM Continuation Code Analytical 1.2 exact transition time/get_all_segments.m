function S =  get_all_segments(flag, xct, dxct, kt, kn, mu, Ft, H)
    Nx = size(flag, 2);
    S = {};
    for i = 1:Nx
        flagT1(:, 1) = flag(1, i, :);
        flagT2(:, 1) = flag(2, i, :);
        xt1 = xct(:, 3 * i - 2);
        xt2 = xct(:, 3 * i - 1);
         xn = xct(:, 3 * i);
        dxt1 = dxct(:, 3 * i - 2);
        dxt2 = dxct(:, 3 * i - 1);
         dxn = dxct(:, 3 * i);
        ft1 = Ft(:, 3 * i - 2);
        ft2 = Ft(:, 3 * i - 1);
        S{1, i} = get_integral_time_position(flagT1, xt1, xn, dxt1, dxn, kt(1, i), kn(i), mu(1, i), ft1, H);
        S{2, i} = get_integral_time_position(flagT2, xt2, xn, dxt2, dxn, kt(2, i), kn(i), mu(2, i), ft2, H);
    end
end