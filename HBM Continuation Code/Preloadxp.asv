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
% * Rx: reaction forces of cpntact part
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
% convert fc to struct type for calling mex function
fc.kn = params.func.fc.kn;
fc.xn0 = params.func.fc.xn0;
fc.mu = params.func.fc.mu;
fc.kt = params.func.fc.kt;
fc.w = params.func.fc.w;

maxiter = params.Newton.maxiter;
epsx = params.Newton.epsx;
epsf = params.Newton.epsf;
Kxx = params.func.CB_MK.Kxx;
% gxpstruct = g(xp', params.func.fc); %%%
gxpstruct = g_mex(xp', fc); 
f = Kxx * xp + gxpstruct.F' + Rx; %%% transpose
w = gxpstruct.w; % update w
disp('preload iterations:'); 
disp('i xp w gxp errxp errf'); 
for i = 1:maxiter
    % Jf = Kxx + finite_diff_jac(@(xp) g(xp, params.func.fc).F, xp'); %%
    Jf = Kxx + finite_diff_jac(@(xp) g_mex(xp, fc).F, xp');
    dxp = - Jf \ f;
    xp = xp + dxp;
    %%%%%
    % gxpstruct = g(xp', params.func.fc); %%%
    gxpstruct = g_mex(xp', fc); 
    f = Kxx * xp + gxpstruct.F' + Rx; %%%
    w = gxpstruct.w; % update w
    %%%%%
    errorxp = norm(dxp);
    errorf = norm(f);
    disp([i, xp(1:3)', w(:,1)', gxpstruct.F(1:3), errorxp, errorf]);
    
    if (errorxp < epsx) && (errorf < epsf)
        gxp = gxpstruct.F';
        return
    end
    if i >= maxiter
        error('Preload displacement of contact part xp is not converged, please check.');
    end
end

end

