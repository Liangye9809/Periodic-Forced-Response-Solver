function [JNL, Mft, JNLt, F, wt] = HBMJACOB_analytical_gf_2dofs_2(xt, kn, xn0, mu, kt, w_in, H, N, nloop)

    [F, wt, Mft] = gf_2dofs(xt, kn, xn0, mu, kt, w_in, nloop);

    
    %%
    a(:, 1) = Mft(1, 1, end-N+1:end);
    b(:, 1) = Mft(2, 1, end-N+1:end);
    c(:, 1) = Mft(3, 1, end-N+1:end);

    JNLt = zeros(2 * N, 4 * H + 2);
    JNLt(1:N, 1) = 0.5 * kt .* c .* a;
    JNLt(1:N, 2 * H + 2) = 0.5 * mu * kn .* c .* b;
    JNLt(N + 1:end, 2 * H + 2) = 0.5 * kn .* c;
    t = ((0:(N-1)) * 2 * pi / N)';
    for i = 1:H
        JNLt(1:N, 2 * i) = kt .* cos(i * t) .* c .* a;
        JNLt(1:N, 2 * i + 1) = kt .* sin(i * t) .* c .* a;

        JNLt(1:N, 2 * H + 1 + 2 * i) = mu * kn .* cos(i * t) .* c .* b;
        JNLt(1:N, 2 * H + 1 + 2 * i + 1) = mu * kn .* sin(i * t) .* c .* b;

        JNLt(N + 1:end, 2 * H + 1 + 2 * i) = kn .* cos(i * t) .* c;
        JNLt(N + 1:end, 2 * H + 1 + 2 * i + 1) = kn .* sin(i * t) .* c;
    end


    JNL = zeros(2 * (2*H+1), 2 * (2*H+1));
    [E, EH] = fft_matrices(N, H);
    % add missing part

    if ismember(0, a)
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
        JNLt(1:N, 1) = JNLt(1:N, 1) - a .* 0.5 .* kt;
        JNLt(1:N, 2 * H + 2) = JNLt(1:N, 2 * H + 2) + c_xn .* 0.5 .* kn .* mu;
        for i = 1:H
            JNLt(1:N, 2 * i) = JNLt(1:N, 2 * i) - a .* kt .* cos(record ./ N .* 2.*pi .* i);
            JNLt(1:N, 2 * i + 1) = JNLt(1:N, 2 * i + 1) - a .* kt .* sin(record ./ N .* 2.*pi .* i);

            JNLt(1:N, 2 * H + 1 + 2 * i) = JNLt(1:N, 2 * H + 1 + 2 * i) + c_xn .* mu .* kn .* cos(record ./ N .* 2.*pi .* i);
            JNLt(1:N, 2 * H + 1 + 2 * i + 1) = JNLt(1:N, 2 * H + 1 + 2 * i + 1) + c_xn .* mu .* kn .* sin(record ./ N .* 2.*pi .* i);
        end

        
    end
    JNL(1:2*H+1, 1:2*H+1) = EH * JNLt(1:N, 1:2*H+1);
    JNL(1:2*H+1, 2*H+2:end) = EH * JNLt(1:N, 2*H+2:end);
    JNL(2*H+2:end, 2*H+2:end) = EH * JNLt(N + 1:end, 2*H+2:end);

end


