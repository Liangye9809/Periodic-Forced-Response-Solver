function JNL = HBMJACOB_analytical(xt, kn, xn0, mu, kt, w_in)

    [N, ~] = size(xt);
    x = xt';
    F = zeros(size(x));

    for j = 1:2
        for i = 1:N
            [F(:, i), w] = gf(x(:, i), kn, xn0, mu, kt, w_in);
            w_in = w; 
        end
    end

    F = F';




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
            C = [0; 1]; % slip
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




