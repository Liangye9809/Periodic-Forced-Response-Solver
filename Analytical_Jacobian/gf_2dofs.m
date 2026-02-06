function [F, wt, Mft] = gf_2dofs(xt, kn, xn0, mu, kt, w_in, nloop)

    [N, M] = size(xt); 
    x = xt;
    F = zeros(size(x));
    wt = zeros(N, 1);
    Mft = zeros(3, 1, N);
    for j = 1:nloop
        for i = 1:N
            [F(i, :), wtemp, Mf] = gf_2dofs_instant(x(i, :), kn, xn0, mu, kt, w_in);
            w_in = wtemp; 
            wt(i) = wtemp;
            Mft(:, :, i) = Mf;
        end
    end

end

