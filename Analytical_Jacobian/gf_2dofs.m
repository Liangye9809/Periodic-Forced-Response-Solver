function [F, wt, Mft] = gf_2dofs(xt, kn, xn0, mu, kt, w_in, nloop)

    [N, M] = size(xt); 
    x = xt;
    F = zeros(nloop * N, 2);
    wt = zeros(nloop * N, 1);
    Mft = zeros(3, 1, N * nloop);
    for j = 1:nloop
        for i = 1:N
            [F((j - 1) * N + i, :), wtemp, Mf] = gf_2dofs_instant(x(i, :), kn, xn0, mu, kt, w_in);
            w_in = wtemp; 
            wt((j - 1) * N + i) = wtemp;
            Mft(:, :, (j - 1) * N + i) = Mf;
        end
    end

end

