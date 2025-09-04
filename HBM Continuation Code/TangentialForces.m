
function [T, w] = TangentialForces(xt, wt, kt, mu, FN)
    if FN > 0
        T = kt * (xt - wt);
        if abs(T) < mu * FN
            w = wt;
        else
            sg = sign(T);
            T = sg * mu * FN;
            w = xt - sg * mu * FN / kt;
        end
    else
        T = zeros(size(xt));
        w = xt;
    end
end

% assume inputs are all row vectors
% function [T, w] = TangentialForces(xt, wt, kt, mu, FN)
%     T_raw = kt .* (xt - wt);
%     slipping = abs(T_raw) >= mu .* FN;
%     sticking = ~slipping;
% 
%     T = T_raw;
%     w = wt;
% 
%     % Stick condition
%     w(sticking) = wt(sticking);
% 
%     % Slip condition
%     sg = sign(T_raw(slipping));
%     T(slipping) = sg .* mu(slipping) .* FN(slipping);
%     w(slipping) = xt(slipping) - T(slipping) ./ kt(slipping);
% end
