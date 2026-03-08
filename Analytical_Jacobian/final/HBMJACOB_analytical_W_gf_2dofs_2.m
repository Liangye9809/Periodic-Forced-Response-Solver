function [JNL, Mft, F, wt] = HBMJACOB_analytical_W_gf_2dofs_2(xt, dx, kn, xn0, mu, kt, w_in, H, N, nloop)

    [F, wt, Mft] = gf_2dofs(xt, kn, xn0, mu, kt, w_in, nloop);

    JNL = zeros(2 * H + 1, 2 * H + 1);
    %%
    a(:, 1) = Mft(1, 1, end-N+1:end);
    b(:, 1) = Mft(2, 1, end-N+1:end);
    c(:, 1) = Mft(3, 1, end-N+1:end);

    



end


