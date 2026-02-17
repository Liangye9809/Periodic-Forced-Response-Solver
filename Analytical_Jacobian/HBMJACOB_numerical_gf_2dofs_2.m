function [JNL, JNLt] = HBMJACOB_numerical_gf_2dofs_2(x, kn, xn0, mu, kt, w_in, H, N, nloop, h, order)
    
    [E, EH] = fft_matrices(N, H);
    JNL = zeros(2 * (2*H+1), 2 * (2*H+1));
    JNLt = finite_diff_jac(@(x) fftFt(x, kn, xn0, mu, kt, w_in, H, N, nloop), x, h, order);
    JNL(1:2*H+1, 1:2*H+1) = EH * JNLt(1:N, 1:2*H+1);
    JNL(1:2*H+1, 2*H+2:end) = EH * JNLt(1:N, 2*H+2:end);
    JNL(2*H+2:end, 2*H+2:end) = EH * JNLt(N + 1:end, 2*H+2:end);

    


end


function [Ft, wt, Mft] = fftFt(x, kn, xn0, mu, kt, w_in, H, N, nloop)

    [E, EH] = fft_matrices(N, H);
    X(:, 1) = x(1:(2 * H + 1));
    X(:, 2) = x(1 + (2 * H + 1):end);
    xt = E * X; 

    [Ft, wt, Mft] = gf_2dofs(xt, kn, xn0, mu, kt, w_in, nloop);
    
    Ft_final = Ft(end - N + 1:end, :);

    Ft = Ft_final(:);


end
