function [JNL, Mft, Gp, dgt, F, wt] = HBMJACOB_analytical_gf_2dofs(xt, kn, xn0, mu, kt, w_in, H, N, nloop)

    [F, wt, Mft] = gf_2dofs(xt, kn, xn0, mu, kt, w_in, nloop);

    
    %%
    [Gp, dgt] = buildGp(Mft(:,:,end-N+1:end), kn, mu, kt, H);

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

        % add missing part

        if ismember(0,Mft(1,1,end-N+1:end))
            [E, EH] = fft_matrices(N, H);
            a(:, 1) = Mft(1,1,end-N+1:end);
            b(:, 1) = Mft(2,1,end-N+1:end);
            c(:, 1) = Mft(3,1,end-N+1:end);
            record = zeros(size(a));
            p = 0;
            keep = 0;
            for i = 1:N

                if i == 1 && a(1) ~=0 
                    for j = 1:N
                        k = N - j + 1;
                        if a(k) == 0
                            p = mod(k, N) + 1;
                            keep = 1;
                            break;
                        end
                    end
                elseif a(i) ~= 0 && keep == 0
                    p = i;
                    keep = 1;
                end

                if a(i) == 0 && keep == 1
                    keep = 0;
                    p = 0;
                end

                record(i) = p;

            end
            c_xn = c .* a .* b(mod(record-2, N));

            dTdXt_time = E * JNL(1:2*H+1, 1:2*H+1);
            dTdXn_time = E * JNL(1:2*H+1, 2*H+2:end);

            dTdXt_time(:, 1) = dTdXt_time(:, 1) - a .* 0.5 .* kt;
            dTdXn_time(:, 1) = dTdXn_time(:, 1) + c_xn .* 0.5 .* kn .* mu;

            for i = 1:H
                dTdXt_time(:, 2 * i) = dTdXt_time(:, 2 * i) - a .* kt .* cos(record ./ N .* 2.*pi .* i);
                dTdXt_time(:, 2 * i + 1) = dTdXt_time(:, 2 * i + 1) - a .* kt .* sin(record ./ N .* 2.*pi .* i);

                dTdXn_time(:, 2 * i) = dTdXn_time(:, 2 * i) + c_xn .* mu .* kn .* cos(record ./ N .* 2.*pi .* i);
                dTdXn_time(:, 2 * i + 1) = dTdXn_time(:, 2 * i + 1) + c_xn .* mu .* kn .* sin(record ./ N .* 2.*pi .* i);
            end

            JNL(1:2*H+1, 1:2*H+1) = EH * dTdXt_time;
            JNL(1:2*H+1, 2*H+2:end) = EH * dTdXn_time;
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





