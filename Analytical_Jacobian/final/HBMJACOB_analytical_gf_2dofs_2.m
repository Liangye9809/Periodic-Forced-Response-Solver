function [JNL, Mft, JNLt, F, wt] = HBMJACOB_analytical_gf_2dofs_2(xt, dx, kn, xn0, mu, kt, w_in, H, N, nloop)

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

    if ismember(0, a) % slip or separation happene

        record_minus = zeros(size(a));
        record_plus = zeros(size(a));
        p_plus = 1; % set default value to 1 not 0, to avoid index of vector irreasonable in code
        p_minus = 1;
        p = 1;
        keep = 0;
        for i = 1:N
        
            if i == 1 && a(1) ~=0 % if begin with stick, means transition time located previously 
                for j = 1:N % backwards loop to find slip to stick transition
                    k = N - j + 1;
                    if a(k) == 0
                        p_plus = mod(k, N + 1) + 1; % t+ transition instant, stick. add last point to create close the cycle [1, N+1]
                        p_minus = mod(k - 1, N) + 1; % t- transition instant, slip
                        p = 0.5 * (p_plus + p_minus); % middle transition instant.
                        keep = 1;
                        break;
                    end
                end
            elseif a(i) ~= 0 && keep == 0
                p_plus = i; % t+ transition instant.
                p_minus = i - 1; % t- transition instant.
                p = 0.5 * (p_plus + p_minus); % middle transition instant.
                keep = 1;
            end
        
            if a(i) == 0 && keep == 1
                keep = 0;
                p_plus = 1;
                p_minus = 1;
                p = 1;
            end
        
            % record(i) = p_plus;
            % record(i) = p;
            record_minus(i) = p_minus;
            record_plus(i) = p_plus;
        
        end
        t_s = (record_minus - 1) ./ N .* 2.*pi; % consider [0, N-1], so there is record - 1 
        c_xn = c .* a .* b(mod(record_plus - 2, N) + 1); % stick part in slip to stick
        dxt = dx(:, 1);
        dxn = dx(:, 2);
        cg_xn = (~c(mod(record_plus - 2, N) + 1)) .* a .* dxt(record_minus) ./ dxn(record_minus); % stick part in gap to stick

        JNLt(1:N, 1) = JNLt(1:N, 1) - a .* 0.5 .* kt;
        
        JNLt(1:N, 2 * H + 2) = JNLt(1:N, 2 * H + 2) + c_xn .* 0.5 .* kn .* mu; % slip to stick correction
        JNLt(1:N, 2 * H + 2) = JNLt(1:N, 2 * H + 2) + kt .* cg_xn .* 0.5; % gap to stick correction

        for i = 1:H
            JNLt(1:N, 2 * i) = JNLt(1:N, 2 * i) - a .* kt .* cos(t_s .* i);
            JNLt(1:N, 2 * i + 1) = JNLt(1:N, 2 * i + 1) - a .* kt .* sin(t_s .* i);

            % slip to stick correction
            JNLt(1:N, 2 * H + 1 + 2 * i) = JNLt(1:N, 2 * H + 1 + 2 * i) + c_xn .* mu .* kn .* cos(t_s .* i);
            JNLt(1:N, 2 * H + 1 + 2 * i + 1) = JNLt(1:N, 2 * H + 1 + 2 * i + 1) + c_xn .* mu .* kn .* sin(t_s .* i);

            % gap to stick correction
            JNLt(1:N, 2 * H + 1 + 2 * i) = JNLt(1:N, 2 * H + 1 + 2 * i)  + kt .* cg_xn .* cos(t_s .* i);
            JNLt(1:N, 2 * H + 1 + 2 * i + 1) = JNLt(1:N, 2 * H + 1 + 2 * i + 1)  + kt .* cg_xn .* sin(t_s .* i);
        end

        
    end

    
    JNL(1:2*H+1, 1:2*H+1) = EH * JNLt(1:N, 1:2*H+1);
    JNL(1:2*H+1, 2*H+2:end) = EH * JNLt(1:N, 2*H+2:end);
    JNL(2*H+2:end, 2*H+2:end) = EH * JNLt(N + 1:end, 2*H+2:end);

end


