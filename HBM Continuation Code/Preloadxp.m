function [xp, gxp] = Preloadxp(Rx, params)
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

% xp = p.xp0;
% gxp = g(xp');
% gxp = gxp';
% f = Kxx * xp + gxp + Rx;
% for i = 1:p.itermax
%     Jf = Kxx + finite_diff_jac(@g,xp');
%     dxp = - Jf \ f;
%     xp = xp + dxp;
%     %%%%%
%     gxp = g(xp');
%     gxp = gxp';
%     f = Kxx * xp + gxp + Rx;
%     %%%%%
%     errorxp = norm(dxp);
%     errorf = norm(f);
%     if (errorxp < p.epsx) && (errorf < p.epsf)
%         return
%     end
% end



xp = params.func.static.xp0;
if size(xp,1) == 1
    xp = xp * ones(size(Rx));
end

maxiter = params.Newton.maxiter;
epsx = params.Newton.epsx;
epsf = params.Newton.epsf;
Kxx = params.func.CB_MK.Kxx;
gxp = g(xp', params.func); %%%
f = Kxx * xp + gxp' + Rx; %%% transpose
disp('preload iterations:'); 
disp('i xp w gxp errxp errf'); 
for i = 1:maxiter
    Jf = Kxx + finite_diff_jac(@(xp) g(xp, params.func), xp'); %%
    dxp = - Jf \ f;
    xp = xp + dxp;
    %%%%%
    gxp = g(xp', params.func); %%%
    f = Kxx * xp + gxp' + Rx; %%%
    %%%%%
    errorxp = norm(dxp);
    errorf = norm(f);
    disp([i,xp(1:3)',params.func.fc.w(:,1)',gxp(1:3),errorxp,errorf]);
    
    if (errorxp < epsx) && (errorf < epsf)
        gxp = gxp';
        return
    end
    if i >= maxiter
        error('Preload displacement of contact part xp is not converged, please check.');
    end
end

end

