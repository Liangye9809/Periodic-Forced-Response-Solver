
function [F, w] = gf(x, fc) % x is a 3*Np dof row vector, represent 2*Np tangential directions and 1*Np normal direction
    
    % fc = fc.fc;
    
    kn = fc.kn; % stiffness in normal direction
    xn0 = fc.xn0; % initial gap between this contact point
    mu = fc.mu; % friciton coefficient of coulomb friction
    kt = fc.kt; % stiffness of 2 tangential direction in column vector
    w = fc.w;
    

    Np = size(kn, 2);
    FN = zeros(1, Np);
    T1 = zeros(1, Np);
    T2 = zeros(1, Np);
    F.F = zeros(3*Np, 1);
    for i = 1:Np
        FN(i) = NormalForces(x(3*(i-1)+3), kn(i), xn0(i));
        [T1(i), w(1, i)] = TangentialForces(x(3*(i-1)+1), w(1, i), kt(1, i), mu(1, i), FN(i));
        [T2(i), w(2, i)] = TangentialForces(x(3*(i-1)+2), w(2, i), kt(2, i), mu(2, i), FN(i));
    end

    FF = [T1; T2; FN];
    F.F = FF(:);
    F.w = w;
end