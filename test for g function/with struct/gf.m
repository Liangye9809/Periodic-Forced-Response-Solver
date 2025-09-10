



function [F, w] = gf(x, fc) % x is a 3*Np dof row vector, represent 2*Np tangential directions and 1*Np normal direction



    kn = fc.kn; % stiffness in normal direction
    xn0 = fc.xn0; % initial gap between this contact point
    mu = fc.mu; % friciton coefficient of coulomb friction
    kt = fc.kt; % stiffness of 2 tangential direction in column vector
    w = fc.w;


    Np = size(kn, 1);
 
    F = zeros(size(x));
    for i = 1:Np
        indx1 = 3 * (i - 1) + 1;
        indx2 = 3 * (i - 1) + 2;
        indx3 = 3 * (i - 1) + 3;
        F(indx3) = NormalForces(x(indx3), kn(i), xn0(i));
        [F(indx1), w(i, 1)] = TangentialForces(x(indx1), w(i, 1), kt(i, 1), mu(i, 1), F(indx3));
        [F(indx2), w(i, 2)] = TangentialForces(x(indx2), w(i, 2), kt(i, 2), mu(i, 2), F(indx3));
    end

end