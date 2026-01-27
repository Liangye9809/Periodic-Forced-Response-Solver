% x are in frequency domain 
% x = [a¹0,a¹1,b¹1,a¹2,b¹2,...,a¹H,b¹H,  a²0,a²1,b²1,a²2,b²2,...,a²H,b²H,...]'

%% original structure
function [F, w] = fftgx_f(x, pfunc) % x(t) = E*X
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
    % [Ft, w] = g_mex(xt + xp', pfunc.fc); % call mex function
    [Ft, w] = g(xt + xp', pfunc.fc); % call mex function inside g
    
    Ft = Ft - gxp';
    hndn = EH * Ft;
    F = hndn(:);
    
end


