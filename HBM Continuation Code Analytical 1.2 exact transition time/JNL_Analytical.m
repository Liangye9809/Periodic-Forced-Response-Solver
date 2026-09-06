function JNL = JNL_Analytical(segments, H, kt, kn, mu) % x is the size of N*3Nx

    Nx = length(segments) / 2; % contact points including T1 and T2
    JNL = zeros(3 * Nx * (2 * H + 1), 3 * Nx * (2 * H + 1));

    for i = 1:Nx
        segmentsT1 = segments{2 * i - 1};
        segmentsT2 = segments{2 * i};
        JNLi = JNL_one_Nx(H, kt(:, i), kn(i), mu(:, i), segmentsT1, segmentsT2);
        indx1 = 3 * (i - 1) * (2 * H + 1) + 1;
        indx2 = 3 * i * (2 * H + 1);
        JNL(indx1:indx2, indx1:indx2) = JNLi;
    end
            

end


function JNLi = JNL_one_Nx(H, kt, kn, mu, segmentsT1, segmentsT2)

    JNLi = zeros(3 * (2 * H + 1), 3 * (2 * H + 1));

    [dF1dX1, dF1dXn, dFndXn] = get_dFdX(segmentsT1, H, kt(1), kn, mu(1));
    
    [dF2dX2, dF2dXn, ~] = get_dFdX(segmentsT2, H, kt(2), kn, mu(2));

    JNLi(1:2 * H + 1, 1:2 * H + 1) = dF1dX1;
    JNLi(1:2 * H + 1, 2 * (2 * H + 1) + 1:end) = dF1dXn;

    JNLi((2 * H + 1) + 1:2 * (2 * H + 1), (2 * H + 1) + 1:2 * (2 * H + 1)) = dF2dX2;
    JNLi((2 * H + 1) + 1:2 * (2 * H + 1), 2 * (2 * H + 1) + 1:end) = dF2dXn;

    JNLi(2 * (2 * H + 1) + 1:end, 2 * (2 * H + 1) + 1:end) = dFndXn;
end


function [dFdX, dFdXn, dFndXn] = get_dFdX(segmentsT, H, kt, kn, mu)
    dFdX   = zeros(2 * H + 1, 2 * H + 1);
    dFdXn  = zeros(2 * H + 1, 2 * H + 1);
    dFndXn = zeros(2 * H + 1, 2 * H + 1);

    k = length(segmentsT);

    if k == 1 % pure stick or whole gap
        if segmentsT(1).value == 2 % pure stick
            dFdX   = kt .* eye(2 * H + 1);
            dFndXn = kn .* eye(2 * H + 1);
            return;
        elseif segmentsT(1).value == 0 % whole gap
            return;
        else
            error('whole flag is %d', segmentsT.value);
        end
    end

    for i = 1:k
        t1 = segmentsT(i).t_start;
        t2 = segmentsT(i).t_end;
        MW = segmentsT(i).MW;
        Mw = segmentsT(i).vw;
        switch segmentsT(i).value
            case 2 % stick after slip or gap
                % [MW, Mw] = fW(t1, t2, H);
                c_vec = c_vector(t1, H);
                dFdX  = dFdX + kt .* MW - Mw * (kt .* c_vec);
                dFndXn = dFndXn + kn .* MW; % contact
                switch segmentsT(mod(i - 2, k) + 1).value

                    case 0 % gap to stick
                        if isempty(segmentsT(i).dxdn)
                            error('no dxdn stored in gap to stick transition!');
                        end
                        dxdn  = segmentsT(i).dxdn;
                        dFdXn = dFdXn + Mw * (kt .* dxdn .* c_vec);

                    case {-1, 1} % slip to stick
                        dFdXn = dFdXn + Mw * (segmentsT(mod(i - 2, k) + 1).value .* mu .* kn .* c_vec);
                end 

            case {-1, 1} % slip
                % [MW, ~] = fW(t1, t2, H);
                dFdXn   = dFdXn + segmentsT(i).value .* mu .* kn .* MW;
                dFndXn  = dFndXn + kn .* MW; % contact
        end

    end

end