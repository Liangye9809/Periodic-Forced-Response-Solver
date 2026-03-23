function [Ft, wt, Mft, dxdnt] = fftFt(x, kn, xn0, mu, kt, w_in, H, N, nloop)

    [E, EH] = fft_matrices(N, H);
    X(:, 1) = x(1:(2 * H + 1));
    X(:, 2) = x(1 + (2 * H + 1):end);
    xt = E * X; 

    [Ft, wt, Mft, dxdnt] = gf_2dofs(xt, kn, xn0, mu, kt, w_in, nloop);
    
    Ft_final = Ft(end - N + 1:end, :);

    Ft = Ft_final(:);


end
