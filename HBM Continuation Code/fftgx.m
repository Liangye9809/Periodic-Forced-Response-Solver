% x are in frequency domain 
% x = [a¹0,a¹1,b¹1,a¹2,b¹2,...,a¹H,b¹H,  a²0,a²1,b²1,a²2,b²2,...,a²H,b²H,...]'
function hndn = fftgx(x, pfunc) % x(t) = E*X
    E = pfunc.HBM.E;
    EH = pfunc.HBM.EH;
    xp = pfunc.static.preload.xp;
    gxp = pfunc.static.preload.gxp;
    n = size(E, 2); % 2H+1
    a = size(x, 1) / n; % number of DOF
    X = zeros(n, a);
    for i = 1:a
        r1 = n * (i - 1) + 1;
        r2 = n * i;
        X(:,i) = x(r1:r2); % reorder in dofs in column
    end
    xt = E * X; 
    % xpt = ones(size(xt)) * diag(xp);
    gt = g(xt + xp', pfunc) - gxp';
    hndn = EH * gt;
    % save Ff_omega1.2.mat hndn
    hndn = hndn(:);
end

