function JNL = HBMJACOB_analytical_gf(xt, kn, xn0, mu, kt, w_in, H)

    [N, M] = size(xt); % N is number of time steps, M is number of dofs - 3Nc
    Nc = M / 3;
    x = xt';
    F = zeros(size(x));
    
    Mft = zeros(5, Nc, N);
    for j = 1:2
        for i = 1:N
            [F(:, i), w, Mf] = gf(x(:, i), kn, xn0, mu, kt, w_in);
            w_in = w; 
            Mft(:, :, i) = Mf;
        end
    end

    % F = F';
    dGmn = buildGp(Mft, kn, mu, kt, H);

    JNL = zeros(3*Nc * (2*H+1), 3*Nc * (2*H+1));

    skip = [1 2; 2 1; 3 1; 3 2];
    for k = 1:Nc

        for i = 1:3
            for j = 1:3
                if ismember([i j], skip, 'rows')
                    continue;
                end
                inda1 = (3 * (k - 1) + i - 1) * (2 * H + 1) + 1;
                inda2 = (3 * (k - 1) + i) * (2 * H + 1);
                indb1 = (3 * (k - 1) + j - 1) * (2 * H + 1) + 1;
                indb2 = (3 * (k - 1) + j) * (2 * H + 1);
                dGtemp(:, 1) = dGmn(3 * (k - 1) + i, 3 * (k - 1) + j, :);
                JNL(inda1:inda2, indb1:indb2) = HBMJACOB_analytical(dGtemp, H);
            end
        end

    end

end

function Gp = buildGp(Mft, kn, mu, kt, H)

    [~, Nc, N] = size(Mft);
    dg = zeros(3*Nc, 3, N);
    [E, EH] = HBM.fft_matrices(N, 2*H);
    for i = 1:N
        for j = 1:Nc
            dg(3 * j - 2, :, i) = Mft(5, j, i) * Mft(1, j, i) * [kt(1, j), 0, 0]...
                                + Mft(5, j, i) * Mft(2, j, i) * [0, 0, mu(1, j) * kn(j)];
            dg(3 * j - 1, :, i) = Mft(5, j, i) * Mft(3, j, i) * [kt(2, j), 0, 0]...
                                + Mft(5, j, i) * Mft(4, j, i) * [0, 0, mu(2, j) * kn(j)];
            dg(3 * j, :, i) = Mft(5, j, i) * [0, 0, kn(j)];
        end
    end
        
    Gp = zeros(3*Nc, 3*Nc, 4*H+1);
    dgj = zeros(N, 1); % store the sigle dg
    Gpj = zeros(4 * H + 1, 1);

    for j = 1:Nc
        % g11
        dgj(:, 1) = dg(3 * j - 2, 1, :);
        Gpj(:, 1) = EH * dgj;
        Gp(3 * j - 2, 3 * j - 2, :) = Gpj;
        % g13
        dgj(:, 1) = dg(3 * j - 2, 3, :);
        Gpj(:, 1) = EH * dgj;
        Gp(3 * j - 2, 3 * j, :) = Gpj;
        % g22
        dgj(:, 1) = dg(3 * j - 1, 2, :);
        Gpj(:, 1) = EH * dgj;
        Gp(3 * j - 1, 3 * j - 1, :) = Gpj;
        % g23
        dgj(:, 1) = dg(3 * j - 1, 2, :);
        Gpj(:, 1) = EH * dgj;
        Gp(3 * j - 1, 3 * j, :) = Gpj;
        % g33
        dgj(:, 1) = dg(3 * j, 3, :);
        Gpj(:, 1) = EH * dgj;
        Gp(3 * j, 3 * j, :) = Gpj;
    end

end

function [F, w, Mf] = gf(x, kn, xn0, mu, kt, w_in) % x is a 3*Np dof row vector, represent 2*Np tangential directions and 1*Np normal direction

    Np = size(kn, 2); % contact number
    F = zeros(size(x));
    Mf = zeros(5, Np); % condition of friction matix
    w = w_in;
    for i = 1:Np
        indx1 = 3 * (i - 1) + 1;
        indx2 = 3 * (i - 1) + 2;
        indx3 = 3 * (i - 1) + 3;
        F(indx3) = NormalForces(x(indx3), kn(i), xn0(i));
        [F(indx1), w(1, i), Mf(1:2, i)] = TangentialForces(x(indx1), w(1, i), kt(1, i), mu(1, i), F(indx3));
        [F(indx2), w(2, i), Mf(3:4, i)] = TangentialForces(x(indx2), w(2, i), kt(2, i), mu(2, i), F(indx3));
        if abs(F(indx3)) > 0
            Mf(5, i) = 1;
        end
    end

end

% C is condition of friction, [stick; slip] 1 is active, 0 is non.
function [T, w, C] = TangentialForces(xt, wt, kt, mu, FN)
    if FN > 0
        T = kt * (xt - wt);
        if abs(T) < mu * FN
            w = wt;
            C = [1; 0]; % stick
        else
            sg = sign(T);
            T = sg * mu * FN;
            w = xt - sg * mu * FN / kt;
            C = sg * [0; 1]; % slip
        end
    else
        T = zeros(size(xt));
        w = xt;
        C = [0; 0]; % gap, no stick nor slip
    end
end


function FN = NormalForces(xn, kn, xn0)
    u = xn - xn0;
    FN = max(0, kn .* u);
end




