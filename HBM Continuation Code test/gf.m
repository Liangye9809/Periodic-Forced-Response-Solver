function [F, w, flag] = gf(x, dx, dx_pre, kn, xn0, mu, kt, w_in, xn_pre, flag_pre, dtheta) % x is a 3*Nx dof row vector, represent 2*Nx tangential directions and 1*Nx normal direction

    Nx = size(kn, 2); % contact number
    F = zeros(size(x));
    w = w_in;
    flag = zeros(2, Nx); % condition of friction matix
    
    for i = 1:Nx
        indx1 = 3 * (i - 1) + 1;
        indx2 = 3 * (i - 1) + 2;
        indx3 = 3 * (i - 1) + 3;
        F(indx3) = NormalForces(x(indx3), kn(i), xn0(i));
        [F(indx1), w(1, i), flag(1, i), ~] = TangentialForces(x(indx1), w(1, i), kt(1, i), mu(1, i), F(indx3), x(indx3), xn_pre(i), kn(i), [dx_pre(indx1), dx(indx1)], [dx_pre(indx3), dx(indx3)], flag_pre(1, i), dtheta);
        [F(indx2), w(2, i), flag(2, i), ~] = TangentialForces(x(indx2), w(2, i), kt(2, i), mu(2, i), F(indx3), x(indx3), xn_pre(i), kn(i), [dx_pre(indx2), dx(indx2)], [dx_pre(indx3), dx(indx3)], flag_pre(2, i), dtheta);
    end

end
