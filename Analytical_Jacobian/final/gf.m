
function [F, w, Mf] = gf(x, kn, xn0, mu, kt, w_in) % x is a 3*Np dof row vector, represent 2*Np tangential directions and 1*Np normal direction

    Nx = size(kn, 2); % contact number
    F = zeros(size(x));
    w = w_in;
    Mf = zeros(5, Nx); % condition of friction matix
    for i = 1:Nx
        indx1 = 3 * (i - 1) + 1;
        indx2 = 3 * (i - 1) + 2;
        indx3 = 3 * (i - 1) + 3;
        F(indx3) = NormalForces(x(indx3), kn(i), xn0(i));
        if F(indx3) > 0
            Mf(5, i) = 1;
        end
        [F(indx1), w(1, i), Mf(1:2, i)] = TangentialForces(x(indx1), w(1, i), kt(1, i), mu(1, i), F(indx3));
        [F(indx2), w(2, i), Mf(3:4, i)] = TangentialForces(x(indx2), w(2, i), kt(2, i), mu(2, i), F(indx3));
    end

end