function [JNL, Mft, Gp, dgt, F, wt] = HBMJACOB_analytical_gf_2dofs(xt, kn, xn0, mu, kt, w_in, H, nloop)

    [F, wt, Mft] = gf_2dofs(xt, kn, xn0, mu, kt, w_in, nloop);

    
    %%
    [Gp, dgt] = buildGp(Mft, kn, mu, kt, H);

    JNL = zeros(2 * (2*H+1), 2 * (2*H+1));

    
        for i = 1:2
            for j = 1:2
                if (i == 2) && (j == 1)
                    continue;
                end
                inda1 = (i - 1) * (2 * H + 1) + 1;
                inda2 = (i) * (2 * H + 1);
                indb1 = (j - 1) * (2 * H + 1) + 1;
                indb2 = (j) * (2 * H + 1);
                dGtemp(:, 1) = Gp(i, j, :);
                JNL(inda1:inda2, indb1:indb2) = HBMJACOB_analytical(dGtemp, H);
            end
        end



end

function [Gp, dgt] = buildGp(Mft, kn, mu, kt, H)

    [~, ~, N] = size(Mft);
    dgt = zeros(2, 2, N);
    [E, EH] = fft_matrices(N, 2*H);
    for i = 1:N
   
        dgt(1, :, i) = Mft(3, 1, i) * Mft(1, 1, i) * [kt, 0]...
                     + Mft(3, 1, i) * Mft(2, 1, i) * [0, mu * kn];
        dgt(2, :, i) = Mft(3, 1, i) * [0, kn];

    end
        
    Gp = zeros(2, 2, 4*H+1);
    dgj = zeros(N, 1); % store the sigle dg
    Gpj = zeros(4 * H + 1, 1);

   
    % g11
    dgj(:, 1) = dgt(1, 1, :);
    Gpj(:, 1) = EH * dgj;
    Gp(1, 1, :) = Gpj;
    % g12
    dgj(:, 1) = dgt(1, 2, :);
    Gpj(:, 1) = EH * dgj;
    Gp(1, 2, :) = Gpj;
    % g22
    dgj(:, 1) = dgt(2, 2, :);
    Gpj(:, 1) = EH * dgj;
    Gp(2, 2, :) = Gpj;


end





