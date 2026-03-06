function [JNL, Mft, F, wt] = HBMJACOB_analytical_W_gf_2dofs_2(xt, dx, kn, xn0, mu, kt, w_in, H, N, nloop)

    [F, wt, Mft] = gf_2dofs(xt, kn, xn0, mu, kt, w_in, nloop);

    JNL = zeros(2 * H + 1, 2 * H + 1);
    %%
    a(:, 1) = Mft(1, 1, end-N+1:end);
    b(:, 1) = Mft(2, 1, end-N+1:end);
    c(:, 1) = Mft(3, 1, end-N+1:end);

    



end


function [C_ss, C_gs, C_s] = get_integral_time_position(a, b, c)
% a is stick condition each time instant, 1 is stick, 0 is non-stick, column vector
% b is slip condition each time instant, 1 or -1 is slip, 0 is non-slip, column vector
% c is contact condition each time instant, 1 is contact, 0 is gap, column vector
% C_ss stick after slip time range
% C_gs stick after gap time range
% C_s slip time range
    N = size(a, 2);
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
end
