function [ft_out, M_fstar] = FFtFactor(ft_in, xt, flag, kt, kn, mu)
    ft_out = ft_in;
    [N, N_3xc] = size(xt);
    M_fstar = zeros(2 * H + 1, N_3xc);
    Nx = N_3xc / 3;
    for i = 1:Nx
        flagT1(:, 1) = flag(1, i, :);
        flagT2(:, 1) = flag(2, i, :);
        xn = xt(:, 3 * i);
        xt1 = xt(:, 3 * i - 2);
        xt2 = xt(:, 3 * i - 1);
        [ft_out(:, 3 * i - 2), M_fstar(:, 3 * i - 2)] = FFtFactor_perT(ft_in(:, 3 * i - 2), xt1, xn, flagT1, kt(1, i), kn(i), mu(1, i));
        [ft_out(:, 3 * i - 1), M_fstar(:, 3 * i - 1)] = FFtFactor_perT(ft_in(:, 3 * i - 1), xt2, xn, flagT2, kt(2, i), kn(i), mu(2, i));
    end
end

function [Ft_out, v_fstar] = FFtFactor_perT(Ft_in, xt, xn, flag, kt, kn, mu) % all the dofs of xt, and Ft
    Ft_out = Ft_in;
    v_fstar = 0;
    N = size(xt, 1);
    % slip and stick transiton
    diffs = [diff(flag); flag(1) - flag(end)];
    ind = ismember(diffs, [1, -1, 3, -3]); % possible slip - stick transition
    trans_ss = find(ind == 1);
    for j = 1:size(trans_ss, 1)
        i_m = trans_ss(j);
        i_p = mod(i_m, N) + 1; % next point
        if flag(i_m) == 2 % stick to slip
            f1_m = Ft_in(i_m); % Ft-
            f1_p = f1_m + kt * (xt(i_p) - xt(i_m));
            f2_p = Ft_in(i_p); % Ft+
            f2_m = sign(f1_m) * mu * kn * xn(i_m);
            [Ft_out(i_m), Ft_out(i_p), f_star] = scaleSlipStick(f1_m, f1_p, f2_m, f2_p);
        end
        % if flag1(i_p) == 2 % slip to stick
        %     f1_m = Ft_in(i_m, 3 * i - 2); % Ft-
        %     f1_p = sign(f1_m) * mu(i) * kn(i) * xn(i_m);
        %     f2_p = Ft_in(i_p, 3 * i - 2); % Ft+
        %     f2_m = f2_p - kt(1, i) * (xt1(i_p) - xt1(i_m));    
        %     [Ft_out(i_m, 3 * i - 2), Ft_out(i_p, 3 * i - 2)] = scaleSlipStick(f1_m, f1_p, f2_m, f2_p);
        % end
        
    end


end

function [Ft_m, Ft_p, f_star] = scaleSlipStick(f1_m, f1_p, f2_m, f2_p)
    fdt = (f2_m - f1_m) / (f1_p - f1_m + f2_m - f2_p);
    f_star = f1_m + (f1_p - f1_m) * fdt;
    % Ft_m = f1_m * 0.5 * (1 + fdt) + 0.5 * (1 - fdt) * f_star;
    % Ft_p = f2_p * 0.5 * (2 - fdt) + 0.5 * fdt * f_star;
    Ft_m = f1_m * 0.5 * (1 + fdt);
    Ft_p = f2_p * 0.5 * (2 - fdt);
end