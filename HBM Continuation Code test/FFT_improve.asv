function dFt = FFT_improve(ft_in, xt, flag, kt, kn, mu)
    dFt = zeros(size(ft_in));
    [N, N_3xc] = size(xt);
    
    Nx = N_3xc / 3;
    for i = 1:Nx
        flagT1(:, 1) = flag(1, i, :);
        flagT2(:, 1) = flag(2, i, :);

        if sum(ismember([-1, 1, 0], flagT1)) > 0 % correction on stick-slip and gap-contact
            xt1 = xt(:, 3 * i - 2);
            ft1 = ft_in(:, 3 * i - 2);
            xn = xt(:, 3 * i);
            dFt1 = FFT_improve_stick_to_slip(ft1, xt1, xn, flagT1, kt(1, i), kn(i), mu(1, i));
            dFt(:, 3 * i - 2) = dFt(:, 3 * i - 2) + dFt1;
        end

        if sum(ismember([-1, 1, 0], flagT2)) > 0 % correction on stick-slip and gap-contact
            xt2 = xt(:, 3 * i - 1);
            ft2 = ft_in(:, 3 * i - 1);
            xn = xt(:, 3 * i);
            dFt2 = FFT_improve_stick_to_slip(ft2, xt2, xn, flagT2, kt(2, i), kn(i), mu(2, i));
            dFt(:, 3 * i - 1) = dFt(:, 3 * i - 1) + dFt2;
        end

        if sum(ismember(0, flagT1)) > 0 && sum(ismember(0, flagT1)) < N % correction on gap-contact
            xn_gap = xt(:, 3 * i);
            ft_gap = ft_in(:, 3 * i - 2:3 * i);
            dFt_gap = FFT_improve_gap_contact(ft_gap, xn_gap);
            dFt(:, 3 * i - 2:3 * i) = dFt(:, 3 * i - 2:3 * i) + dFt_gap;
        end

    end
end


function dFt = FFT_improve_stick_to_slip(ft, xt, xn, flag, kt, kn, mu)
    dFt = zeros(size(ft));
    N = size(xt, 1);
    diffs = [diff(flag); flag(1) - flag(end)];
    ind = ismember(diffs, [1, -1, 3, -3, -2]); % possible stick-slip or stick-gap transition
    trans_ss = find(ind == 1); % transition position
    for j = 1:size(trans_ss, 1)
        i_m = trans_ss(j);
        i_p = mod(i_m, N) + 1; % next point
          if flag(i_m) == 2 % stick to slip or gap
            f1_m = ft(i_m); % Ft-
            f1_p = f1_m + kt * (xt(i_p) - xt(i_m));
            if flag(i_p) == 0 %  stick to gap, possible stick-slip-gap in one interval
                fdt_gap = xn(i_m) / (xn(i_m) - xn(i_p));
                f_gap = f1_m + (f1_p - f1_m) * fdt_gap;
                if f_gap ~= 0 % nonzero force at lift-off implies hidden slip
                    f2_m = sign(f_gap) * mu * kn * xn(i_m);
                    fdt = fdt_gap * (f2_m - f1_m) / (f_gap - f1_m + f2_m);
                    f_star = f1_m + (f1_p - f1_m) * fdt;
                    dFt(i_m) = 0.5 * (fdt - fdt_gap) * f1_m ...
                        + 0.5 * fdt_gap * (1 - fdt) * f_star;
                    dFt(i_p) = 0.5 * fdt_gap * fdt * f_star;
                end
            else % stick to slip
                f2_p = ft(i_p); % Ft+
                f2_m = sign(f1_p) * mu * kn * xn(i_m);
                fdt = (f2_m - f1_m) / (f2_m - f1_m + f1_p - f2_p);
                f_star = f1_m + (f1_p - f1_m) * fdt;
                dFt(i_m) = 0.5 * (1 - fdt) * (f_star - f1_m);
                dFt(i_p) = 0.5 * fdt * (f_star - f2_p);
            end
          end
    end
end

function dFt_gap = FFT_improve_gap_contact(ft, xn)
    dFt_gap = zeros(size(ft));
    N = size(xn, 1);
    s = sign(xn);
    diffs = [diff(s); s(1) - s(end)];

    trans_contact = find(diffs == -2); % contact to gap
    for j = 1:size(trans_contact, 1)
        % for Fn
        i_m = trans_contact(j);
        i_p = mod(i_m, N) + 1; % next point
        fdt = xn(i_m) / (xn(i_m) - xn(i_p));
        dFt_gap(i_m, :) = 0.5 * (fdt - 1) * ft(i_m, :);

    end

    trans_gap = find(diffs == 2); % gap to contact
    for j = 1:size(trans_gap, 1)
        % For Fn
        i_m = trans_gap(j);
        i_p = mod(i_m, N) + 1; % next point
        fdt = xn(i_m) / (xn(i_m) - xn(i_p));
        dFt_gap(i_p, :) = -0.5 * fdt * ft(i_p, :);


    end
end