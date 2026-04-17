function [xp, gxp, w] = Preloadxp(Rx, params)
% *************************************************************************************************
% this function calculate the contact static displacement in preload case 
% by the following nonlinear system using Newton Method
% 
% Kxx * xp + g(xp) + Rx = 0;
% 
% *************************************************************************************************
%   INPUTS:
% * Kxx: Nc*Nc matrix in sparse form
% * Rx: reaction forces of contact part
% * p.xp0: initial displacement in Newton Method
% * p.epsx: xp tolerance in Newton Method
% * p.epsf: f(xp) tolerance in Newton Method
% * p.itermax: maximum iteration steps
% 
%   OUTPUTS:
% * xp: the contact static displacement in preload case
% 
% Written by Liu Liangye on April 07, 2025
% *************************************************************************************************



xp = params.func.static.xp0;
if size(xp,1) == 1
    xp = xp * ones(size(Rx));
end
fc = params.func.fc;

maxiter = params.Newton.maxiter;
epsx = params.Newton.epsx;
epsf = params.Newton.epsf;
Kxx = params.func.CB_MK.Kxx;
% [F, w] = g_mex(xp', fc); 
[F, w] = g(xp', fc); 
f = Kxx * xp + F' + Rx; %%% transpose
disp('preload iterations:'); 
disp('i xp w gxp errxp errf'); 
for i = 1:maxiter
    Jf = Kxx + finite_diff_jac(@(xp) Getgx(xp, fc), xp');
    dxp = - Jf \ f;
    xp = xp + dxp;
    %%%%%
    % [F, w] = g_mex(xp', fc); 
    [F, w] = g(xp', fc); 
    f = Kxx * xp + F' + Rx; %%% transpose
    %%%%%
    errorxp = norm(dxp);
    errorf = norm(f);
    disp([i, xp(1:3)', w(:,1)', F(1:3), errorxp, errorf]);
    
    if (errorxp < epsx) && (errorf < epsf)
        gxp = F';
        return
    end
    if i >= maxiter
        error('Preload displacement of contact part xp is not converged, please check.');
    end
end

end

