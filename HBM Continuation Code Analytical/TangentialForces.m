% function [T, w, Flag] = TangentialForces(xt, wt, kt, mu, FN, xn, xn_pre)
%     if FN > 0
%         T = kt * (xt - wt);
%         if abs(T) < mu * FN
%             w = wt;
%             Flag = 2; % stick
%             if xn_pre <= 0 % previous is gap, means gap to stick
%                 dxdti = (xt - wt) / ((xn - xn_pre) + 1e-16);
%                 w = wt - xn_pre * dxdti; % update w
%                 T = kt * (xt - w); % Does it need to be compare with mu*Fn?
%             end
%         else
%             sg = sign(T);
%             T = sg * mu * FN;
%             w = xt - sg * mu * FN / kt;
%             Flag = sg; % slip
%         end
%     else
%         T = zeros(size(xt));
%         w = xt;
%         Flag = 0; % gap
%     end
% end


function [T, w, Flag] = TangentialForces(xt, wt, kt, mu, FN, xn, xn_pre)
    wp = wt; % previous w
    if FN > 0 % contact
        T = kt * (xt - wp);
        if xn_pre <= 0 % previous is gap, means gap to stick
            dxdti = (xt - wt) / ((xn - xn_pre) + 1e-16);
            wp = wt - xn_pre * dxdti; % update w
            T = kt * (xt - wp);
        end
        if abs(T) < mu * FN
            w = wp;
            Flag = 2; % stick
        else
            sg = sign(T);
            T = sg * mu * FN;
            w = xt - sg * mu * FN / kt;
            Flag = sg; % slip
        end
    else
        T = zeros(size(xt));
        w = xt;
        Flag = 0; % gap
    end
end