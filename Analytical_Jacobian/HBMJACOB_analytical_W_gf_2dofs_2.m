function JNL = HBMJACOB_analytical_W_gf_2dofs_2(dx, kn, mu, kt, H, N, Mft, dxdn)

    % [F, wt, Mft] = gf_2dofs(xt, kn, xn0, mu, kt, w_in, nloop);

    JNL = zeros(2 * (2*H+1), 2 * (2*H+1));
    %%
    a(:, 1) = Mft(1, 1, end-N+1:end);
    b(:, 1) = Mft(2, 1, end-N+1:end);
    c(:, 1) = Mft(3, 1, end-N+1:end);
    [C_ps, C_ss, C_gs, C_sp, C_sm, C_c] = get_integral_time_position(a, b, c);
    
    % pure stick case
    if C_ps == 1
        I = eye(2 * H + 1);
        JNL(1:2 * H + 1, 1:2 * H + 1) = kt .* I;
        JNL(2 * H + 1 + 1:end, 2 * H + 1 + 1:end) = kn .* I;
        return;
    end

    % slip to stick part
    i = 1;
    while C_ss(i, 1) >= 0
        t_ss = C_ss ./ N .* (2 * pi);
        [MW, Mw] = fW1(t_ss(i, 1), t_ss(i, 2), H);
        c_vec = c_vector(t_ss(i, 1), H);

        JNL(1:2 * H + 1, 1:2 * H + 1) = JNL(1:2 * H + 1, 1:2 * H + 1) + kt .* MW - Mw * (kt .* c_vec);
        JNL(1:2 * H + 1, 2 * H + 1 + 1:end) = JNL(1:2 * H + 1, 2 * H + 1 + 1:end) + Mw * (b(mod(ceil(C_ss(i, 1)) - 1, N) + 1) .* mu .* kn .* c_vec);
        
        i = i + 1;
    end

    % gap to stick part
    i = 1;
    while C_gs(i, 1) >= 0
        t_gs = C_gs ./ N .* (2 * pi);
        [MW, Mw] = fW1(t_gs(i, 1), t_gs(i, 2), H);
        c_vec = c_vector(t_gs(i, 1), H);
        JNL(1:2 * H + 1, 1:2 * H + 1) = JNL(1:2 * H + 1, 1:2 * H + 1) + kt .* MW - Mw * (kt .* c_vec);
        
        JNL(1:2 * H + 1, 2 * H + 1 + 1:end) = JNL(1:2 * H + 1, 2 * H + 1 + 1:end) + Mw * (kt .* dxdn(i) .* c_vec);
        i = i + 1;
    end

    % slip minus part
    i = 1;
    while C_sm(i, i) >= 0
        t_sm = C_sm ./ N .* (2 * pi);
        [MW, ~] = fW1(t_sm(i, 1), t_sm(i, 2), H);
        JNL(1:2 * H + 1, 2 * H + 1 + 1:end) = JNL(1:2 * H + 1, 2 * H + 1 + 1:end) - mu .* kn .* MW;
        i = i + 1;
    end

    % slip plus part
    i = 1;
    while C_sp(i, i) >= 0
        t_sp = C_sp ./ N .* (2 * pi);
        [MW, ~] = fW1(t_sp(i, 1), t_sp(i, 2), H);
        JNL(1:2 * H + 1, 2 * H + 1 + 1:end) = JNL(1:2 * H + 1, 2 * H + 1 + 1:end) + mu .* kn .* MW;
        i = i + 1;
    end

    % contact part
    i = 1;
    while C_c(i, i) >= 0
        t_c = C_c ./ N .* (2 * pi);
        [MW, ~] = fW1(t_c(i, 1), t_c(i, 2), H);
        JNL(2 * H + 1 + 1:end, 2 * H + 1 + 1:end) = JNL(2 * H + 1 + 1:end, 2 * H + 1 + 1:end) + kn .* MW;
        i = i + 1;
    end

end


